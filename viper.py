#! /usr/bin/env python3

# ./viper.py data/TLS/betgem/BETA_GEM.fits data/TLS/betgem/pepsib.20150409.000.sxt.awl.all6
# ./viper.py data/TLS/hd189733/TV00001.fits data/TLS/Deconv/HD189733.model
# ./viper.py "data/TLS/hd189733/*" data/TLS/Deconv/HARPS.2006-09-08T02\:12\:38.604_s1d_A.fits

import argparse
import glob
import importlib
import os
import time

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from astropy.io import fits

from gplot import *
gplot.tmp = '$'
gplot.colors('classic')
from pause import pause

#from inst.inst_TLS import Spectrum, Tpl, FTS
#from inst.inst_CES import Spectrum, Tpl, FTS
from model import model, IPs, show_model
import vpr

viperdir = os.path.dirname(os.path.realpath(__file__)) + os.sep

c = 3e5   # [km/s] speed of light

o = 18; lmin = 5240; lmax = 5390

dirname = r''

ftsname = dirname + 'lib/TLS/FTS/TLS_I2_FTS.fits'
obsname = dirname + 'data/TLS/betgem/BETA_GEM.fits'
tplname = dirname + 'data/TLS/betgem/pepsib.20150409.000.sxt.awl.all6'
obsname = dirname + 'data/TLS/hd189733/TV00001.fits'
tplname = dirname + 'data/TLS/Deconv/HD189733.model'
tplname = dirname + 'data/TLS/Deconv/HARPS.2006-09-08T02:12:38.604_s1d_A.fits'

nset = None

def arg2slice(arg):
   """Convert string argument to a slice."""
   # We want four cases for indexing: None, int, list of ints, slices.
   # Use [] as default, so 'in' can be used.
   if isinstance(arg, str):
       arg = eval('np.s_['+arg+']')
   return [arg] if isinstance(arg, int) else arg

def arg2range(arg):
    return  eval('np.r_['+arg+']')

if __name__ == "__main__":
    # check first the instrument with preparsing
    insts = ['TLS', 'CES']
    preparser = argparse.ArgumentParser(add_help=False)
    preparser.add_argument('args', nargs='*')
    preparser.add_argument('-inst', help='Instrument', default='TLS', choices=insts)
    preargs, _ =  preparser.parse_known_args()

    inst = preargs.inst
    inst = importlib.import_module('inst.inst_'+inst)
    FTS = inst.FTS
    Spectrum = inst.Spectrum
    Tpl = inst.Tpl

    parser = argparse.ArgumentParser(description='VIPER - velocity and IP Estimator', add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    argopt = parser.add_argument   # function short cut
    argopt('obspath', help='Filename of observation', default='data/TLS/betgem/BETA_GEM.fits', type=str)
    argopt('tplname', help='Filename of template', default='data/TLS/betgem/pepsib.20150409.000.sxt.awl.all6', type=str)
    argopt('-inst', help='Instrument', default='TLS', choices=insts)
    argopt('-fts', help='Filename of template', default=viperdir + FTS.__defaults__[0], dest='ftsname', type=str)
    argopt('-ip', help='IP model (g: Gaussian, sg: Supergaussian, mg: Multiplegaussian)', default='g', choices=['g', 'sg','mg'], type=str)
    argopt('-look', nargs='?', help='See final fit of chunk', default=[], const=':100', type=arg2range)
    argopt('-lookguess', nargs='?', help='Show inital model', default=[], const=':100', type=arg2range)
    argopt('-lookpar', nargs='?', help='See parameter of chunk', default=[], const=':100', type=arg2range)
    argopt('-nset', help='index for spectrum', default=':', type=arg2slice)
    argopt('-nexcl', help='Pattern ignore', default=[], type=arg2range)
    argopt('-oset', help='index for order', default='18:30', type=arg2slice)
    argopt('-tag', help='Output tag for filename', default='tmp', type=str)
    argopt('-vg', help='RV guess', default=1., type=float)   # slightly offsetted
    argopt('-?', '-h', '-help', '--help',  help='show this help message and exit', action='help')

    args = parser.parse_args()
    globals().update(vars(args))

####  FTS  ####

# using the supersampled log(wavelength) space with knot index j

w_I2, f_I2, xj_full, iod_j_full = FTS(ftsname)

orders = np.r_[oset] # np.arange(18,30)
print(orders)

rv = np.nan * orders
e_rv = np.nan * orders


def fit_chunk(o, obsname):
    ####  data TLS  ####
    w_i, f_i, bp, bjd, berv = Spectrum(obsname, o=o)
    i = np.arange(f_i.size)
    #i_ok = slice(400,1700) # probably the wavelength solution of the template is bad
    mskatm = interp1d(*np.genfromtxt(viperdir+'lib/mask_vis1.0.dat').T)
    bp[mskatm(w_i) > 0.1] |= 16
    i_ok, = np.where(bp==0)

    ####  stellar template  ####
    w_tpl, f_tpl = Tpl(tplname, o=o)

    lmin = max(w_tpl[0], w_i[0], w_I2[0])
    lmax = min(w_tpl[-1], w_i[-1], w_I2[-1])

    # display
    # pre-look raw input
    s = slice(*np.searchsorted(w_I2, [lmin, lmax]))
    s_s = slice(*np.searchsorted(w_tpl, [lmin, lmax]))
    vcut = 30
    sj = slice(*np.searchsorted(xj_full, [np.log(lmin)+vcut/c, np.log(lmax)-vcut/c])) # reduce range by 100 km/s

    # prepare input; convert discrete data to model

    # using the supersampled log(wavelength) space with knot index j
    xj = xj_full[sj]
    iod_j = iod_j_full[sj]

    #pause()
    if 0:
        # plot data, template, and iodine without any modifications
        gplot(w_I2[s], f_I2[s], 'w l lc 9,', w_tpl[s_s], f_tpl[s_s], 'w l lc 3,', w_i, f_i, 'w lp lc 1 pt 7 ps 0.5')
        # plot with some scaling
        #gplot(w_I2[s], f_I2[s]/1.18, 'w l lc 9,', w_tpl[s_s]*(1+12/c), f_tpl[s_s], 'w l lc 3,', w_i, f_i/1.04, 'w lp lc 1 pt 7 ps 0.5')


    # convert discrete template into a function
    S_star = interp1d(np.log(w_tpl)-berv/c, f_tpl)


    # setup the model
    IP = IPs[ip]
    S_mod = model(S_star, xj, iod_j, IP)

    if 0:
        # plot the IP
        gplot(S_mod.vk, S_mod.IP(S_mod.vk))

    if 0:
       # plot again, now the stellar template can be interpolated
       gplot(np.exp(xj), iod_j, S_star(xj), 'w l lc 9, "" us 1:3 w l lc 3')


    # Now wavelength solution

    # mapping between pixel and wavelength

    #lam(x) = b0 + b1 * x + b2 * x^2
    lam = np.poly1d([w_i[0], (w_i[-1]-w_i[0])/w_i.size][::-1])

    # trim the observation to a range valid for the model
    i_ok = slice(*np.searchsorted(np.log(w_i), [xj[0]+100/c, xj[-1]-100/c]))

    i_ok = np.r_[i_ok][bp[i_ok]==0]

    # an initial parameter set
    v = vg   # a good guess for the stellar RV is needed
    a = ag = [np.mean(f_i) / np.mean(S_star(np.log(w_i[i_ok]))) / np.mean(iod_j)] 
    b = bg = [w_i[0], (w_i[-1]-w_i[0])/w_i.size] # [6128.8833940969, 0.05453566108124]
    b = bg = np.polyfit(i[i_ok], w_i[i_ok], 3)[::-1]
    #show_model(i[i_ok], f_i[i_ok], S_b(i[i_ok],*bg), res=False)
    s = sg = [1.] if ip=='g' else [1., 2.]

    if 0:
        # a simple call to the forward model
        # Si_mod = S_mod(i[i_ok], v=0, a=a, b=b, s=s)
        # show the start guess
        S_mod.show([v,a,b,s], i[i_ok], f_i[i_ok], res=False, dx=0.1)
        #pause()

    if 0:
        # A wrapper to fit the continuum
        S_a = lambda x, a0: S_mod(x, v, [a0], b, s)
        a, e_a = curve_fit(S_a, i[i_ok], f_i[i_ok])
        #show_model(i[i_ok], f_i[i_ok], S_a(i[i_ok],*a), res=False)

    if 0:
        # A wrapper to fit the wavelength solution
        S_b = lambda x, b0,b1,b2,b3: S_mod(x, v, a, [b0,b1,b2,b3], s)
        b, e_b = curve_fit(S_b, i[i_ok], f_i[i_ok], p0=bg)
        #show_model(i[i_ok], f_i[i_ok], S_b(i[i_ok], *bg))
        #show_model(i[i_ok], f_i[i_ok], S_b(i[i_ok], *b))
        #gplot+(i[i_ok], S_star(np.log(np.poly1d(b[::-1])(i[i_ok]))+(v)/c), 'w lp ps 0.5')



    if 0:
        # fit v, a and b simulatenously
        S_vab = lambda x, v, a, b0,b1,b2,b3: S_mod(x, v, [a], [b0,b1,b2,b3], 2.2)
        p_vab, e_p = curve_fit(S_vab, i[i_ok], f_i[i_ok], p0=[v, 1, *bg])
        show_model(i[i_ok], f_i[i_ok], S_vab(i[i_ok], *p_vab))

    if 0:
        # fit v, a and b simulatenously
        S_vabs = lambda x, v, a, b0,b1,b2,b3, s: S_mod(x, v, [a], [b0,b1,b2,b3], s)
        p_vabs, e_p_vabs = curve_fit(S_vabs, i[i_ok], f_i[i_ok], p0=[*p_va, 2.2], epsfcn=1e-12)
        #rvo, e_rvo = p_vabs[0], np.diag(e_p_vabs)[0]**0.5
        #print(o, rvo, e_rvo)
        show_model(i[i_ok], f_i[i_ok], S_vabs(i[i_ok], *p_vabs))

    if 'ip' == 'sg':
        S = lambda x, v, a0,a1,a2,a3, b0,b1,b2,b3, s,e: S_mod(x, v, [a0,a1,a2,a3], [b0,b1,b2,b3], [s,e])
    else:
        S = lambda x, v, a0,a1,a2,a3, b0,b1,b2,b3, s: S_mod(x, v, [a0,a1,a2,a3], [b0,b1,b2,b3], [s])

    if o in lookguess:
        pg = [v, a+[0]*3, bg, s]
        prms = S_mod.show(pg, i[i_ok], f_i[i_ok], dx=0.1)
        pause()

    p, e_p = curve_fit(S, i[i_ok], f_i[i_ok], p0=[v]+a+[0]*3+[*bg]+s, epsfcn=1e-12)
    rvo, e_rvo = 1000*p[0], 1000*np.diag(e_p)[0]**0.5   # convert to m/s
    #show_model(i[i_ok], f_i[i_ok], S(i[i_ok], *p_vabs))
    prms = S_mod.show([p[0], p[1:5], p[5:9], p[9:]], i[i_ok], f_i[i_ok], dx=0.1)
    # gplot+(w_tpl[s_s]*(1-berv/c), f_tpl[s_s]*ag, 'w lp lc 4 ps 0.5')

    # error estimation
    # uncertainty in continuum
    xl = np.log(np.poly1d(p[5:9][::-1])(i))
    Cg = np.poly1d(ag[::-1])(xl)        # continuum guess
    Cp = np.poly1d(p[1:5][::-1])(xl)    # best continuum
    X = np.vander(xl,4)[:,::-1].T
    e_Cp = np.einsum('ji,jk,ki->i', X, e_p[1:5,1:5], X)**0.5
    # uncertainty in wavelength solution
    lam_g = np.poly1d(bg[::-1])(i)
    lam_p = np.poly1d(p[5:9][::-1])(i)
    e_lam = np.einsum('ji,jk,ki->i', X, e_p[5:9,5:9], X)**0.5
    e_wavesol = np.sum((e_lam/lam_p*3e8)**-2)**-0.5

    if o in look:
        pause('look ', o)  # globals().update(locals())

    #show_model(i[i_ok], f_i[i_ok], S_b(i[i_ok], *bg))
    #show_model(i[i_ok], f_i[i_ok], S_vabs(i[i_ok], *p))
    #gplot+(i[i_ok], S_star(np.log(np.poly1d(b[::-1])(i[i_ok]))+(v)/c), 'w lp ps 0.5')
    if o in lookpar:
        # compare the wavelength solutions
        #show_model(i, np.poly1d(b[::-1])(i), np.poly1d(bg[::-1])(i), res=True)
        gplot.reset()
        gplot.multiplot("layout 2,2")
        gplot.xlabel('"pixel"').ylabel('"k(x2)"')
        gplot.mxtics().mytics()
        gplot('[][0:]', i, Cg, Cp, e_Cp, 'w l lc 9 t "guess",  "" us 1:3 w l lc 3, "" us 1:($3-$4):($3+$4) w filledcurves fill fs transparent solid 0.2 lc 3 t "1{/Symbol s}" ')
        gplot.xlabel('"pixel"').ylabel('"deviation c * ({/Symbol l} / {/Symbol l}_{guess} - 1) [km/s]"')
        gplot(i, (lam_p/lam_g-1)*c, ((lam_p-e_lam)/lam_g-1)*c, ((lam_p+e_lam)/lam_g-1)*c, 'w l lc 3, "" us 1:3:4 w filledcurves fill fs transparent solid 0.2 lc 3 t "1{/Symbol s}"')
        gplot.xlabel('"[km/s]"').ylabel('"contribution"')
        e_s = e_p[9,9]**0.5
        gplot(S_mod.vk, S_mod.IP(S_mod.vk, *s), ' lc 9 ps 0.5 t "IP_{guess}", ', S_mod.vk, S_mod.IP(S_mod.vk,*p[9:]),  S_mod.IP(S_mod.vk,[[p[9]-e_s],*p[10:]]),  S_mod.IP(S_mod.vk, [p[9]+e_s, *p[10:]]), 'lc 3 ps 0.5 t "IP", "" us 1:3:4 w filledcurves fill fs transparent solid 0.2 lc 3 t "1{/Symbol s}"')
        gplot.unset('multiplot')
        pause()
      
    return rvo, e_rvo, bjd, berv, p, e_p, prms


rvounit = open(tag+'.rvo.dat', 'w')
parunit = open(tag+'.par.dat', 'w')
# file header
print('BJD', 'RV', 'e_RV', 'BERV', *sum(zip(map("rv{}".format, orders), map("e_rv{}".format, orders)),()), 'filename', file=rvounit)
print('BJD', 'o', *sum(zip(map("p{}".format, range(10)), map("e_p{}".format, range(10))),()), 'prms', file=parunit)

obsnames = sorted(glob.glob(obspath))[nset]
obsnames = [x for i,x in enumerate(obsnames) if i not in nexcl]
N = len(obsnames)
if not N: pause('no files: ', obspath)
T = time.time()


for n,obsname in enumerate(obsnames):
    filename = os.path.basename(obsname)
    print("%2d/%d"% (n+1,N), obsname)
    for i_o, o in enumerate(orders):
        gplot.key('title "%s (n=%s, o=%s)"'% (filename, n+1, o))
        rv[i_o], e_rv[i_o], bjd,berv, p, e_p, prms = fit_chunk(o, obsname=obsname)
#        try:
#            rv[i_o], e_rv[i_o], bjd,berv, p, e_p  = fit_chunk(o, obsname=obsname)
#        except Exception as e:
#            if repr(e) == 'BdbQuit()':
#               exit()
#            print("Order failed due to:", repr(e))

        print(n+1, o, rv[i_o], e_rv[i_o])
        print(bjd, o, *sum(zip(p, np.diag(e_p)),()), prms, file=parunit)
        #pause()

    oo = np.isfinite(e_rv)
    if oo.sum() == 1:
        RV = rv[oo][0]
        e_RV = e_rv[oo][0]
    else:
        RV = np.mean(rv[oo])
        e_RV = np.std(rv[oo])/(oo.sum()-1)**0.5
    print('RV:', RV,e_RV, bjd, berv)

    print(bjd, RV, e_RV, berv, *sum(zip(rv, e_rv),()), filename, file=rvounit)
    print(file=parunit)
    #vpr.plot_rvo(rv, e_rv)

rvounit.close()

T = time.time() - T
Tfmt = lambda t: time.strftime("%Hh%Mm%Ss", time.gmtime(t))
print("processing time total:       ", Tfmt(T))
print("processing time per spectrum:", Tfmt(T/N))
print("processing time per chunk:   ", Tfmt(T/N/orders.size))

vpr.plot_RV(tag+'.rvo.dat')
pause()

print('Done.')
