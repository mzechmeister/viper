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
from scipy.optimize import curve_fit
from astropy.io import fits
import astropy.units as u

from gplot import *
gplot.colors('classic')
gplot2 = Gplot()
from pause import pause

from model import model, model_bnd, IPs, show_model
from targ import Targ
import vpr

viperdir = os.path.dirname(os.path.realpath(__file__)) + os.sep

c = 3e5   # [km/s] speed of light

targ = None
modset = {}   # model setting parameters
insts = ['TLS', 'CES', 'KECK', 'UVES', 'OES']

def arg2slice(arg):
    """Convert string argument to a slice."""
    # We want four cases for indexing: None, int, list of ints, slices.
    # Use [] as default, so 'in' can be used.
    if isinstance(arg, str):
        arg = eval('np.s_['+arg+']')
    return [arg] if isinstance(arg, int) else arg

def arg2range(arg):
    return  eval('np.r_['+arg+']')


def SSRstat(vgrid, SSR, dk=1, plot='maybe', N=None):
    '''   
    Analyse chi2 peak.

    Parameters
    ----------
    N: Number of data points in the fit. Needed to estimate 1 sigma uncertainty and when SSR are not chi2 values.

    '''
    k = np.argmin(SSR[dk:-dk]) + dk   # best point (exclude borders)
    vpeak = vgrid[k-dk:k+dk+1]
    SSRpeak = SSR[k-dk:k+dk+1] - SSR[k]
    v_step = vgrid[1] - vgrid[0]
    # interpolating parabola a0+a1*x+a2*x**2 (direct solution) through the three pixels in the minimum
    a = np.array([0, (SSR[k+dk]-SSR[k-dk])/(2*v_step), (SSR[k+dk]-2*SSR[k]+SSR[k-dk])/(2*v_step**2)])  # interpolating parabola for even grid
    v = (SSR[k+dk]-SSR[k-dk]) / (SSR[k+dk]-2*SSR[k]+SSR[k-dk]) * 0.5 * v_step

    v = vgrid[k] - a[1]/2./a[2]   # position of parabola minimum
    e_v = np.nan
    if -1 in SSR:
        print('opti warning: bad ccf.')
    elif a[2] <= 0:
        print('opti warning: a[2]=%f<=0.' % a[2])
    elif not vgrid[0] <= v <= vgrid[-1]:
        print('opti warning: v not in [va,vb].')
    else:
        e_v = 1. / a[2]**0.5
        if N:
            # Rescale the variances (parabola) such that the minimum has value N, implying chi^2_red == 1.
            # Then derive one sigma uncertainty from Delta chi^2 = 1.
            SSRmin = SSR[k] + a[0] - a[1] * a[1]/2./a[2] + a[2] * (a[1]/2./a[2]) **2
            e_v *= (SSRmin/N) **0.5
 

    if (plot==1 and np.isnan(e_v)) or plot==2:
        gplot2.yrange('[*:%f]' % np.max(SSR))
        gplot2(vgrid, SSR-SSR[k], " w lp, v1="+str(vgrid[k])+", %f+(x-v1)*%f+(x-v1)**2*%f," % tuple(a), [v,v], [0,SSR[1]], 'w l t "%f km/s"'%v)
        gplot2+(vpeak, SSRpeak, ' lt 1 pt 6; set yrange [*:*]')
        pause(v)
    return v, e_v, a


if __name__ == "__main__":
    # check first the instrument with preparsing
    preparser = argparse.ArgumentParser(add_help=False)
    preparser.add_argument('args', nargs='*')
    preparser.add_argument('-inst', help='Instrument', default='TLS', choices=insts)
    preargs, _ =  preparser.parse_known_args()

    Inst = importlib.import_module('inst.inst_'+preargs.inst)
    FTS = Inst.FTS
    Tpl = Inst.Tpl
    Spectrum = Inst.Spectrum
    iset = getattr(Inst, 'iset', slice(None))
    oset = getattr(Inst, 'oset')

    parser = argparse.ArgumentParser(description='VIPER - velocity and IP Estimator', add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    argopt = parser.add_argument   # function short cut
    argopt('obspath', help='Filename of observation', default='data/TLS/betgem/BETA_GEM.fits', type=str)
    argopt('tplname', help='Filename of template', nargs='?', type=str)
    argopt('-inst', help='Instrument', default='TLS', choices=insts)
    argopt('-fts', help='Filename of template', default=viperdir + FTS.__defaults__[0], dest='ftsname', type=str)
    argopt('-ip', help='IP model (g: Gaussian, ag: asymmetric (skewed) Gaussian, sg: super Gaussian, mg: multiple Gaussians, bnd: bandmatrix)', default='g', choices=['g', 'ag', 'sg', 'mg', 'bnd'], type=str)
    argopt('-chunks', nargs='?', help='Divide one order into a number of chunks', default=1, type=int)
    argopt('-dega', nargs='?', help='Polynomial degree for flux normalisation.', default=3, type=int)
    argopt('-degb', nargs='?', help='Polynomial degree for wavelength scale l(x).', default=3, type=int)
    argopt('-degc', nargs='?', help='Number of additional parameters.', default=0, const=1, type=int)
    argopt('-demo', nargs='?', help='Demo plots. Use -8 to skip plots 1,2,4).', default=0, const=-1, type=int)
    argopt('-iset', help='maximum', default=iset, type=arg2slice)
    argopt('-kapsig', help='Kappa sigma clipping value', default=None, type=float)
    argopt('-look', nargs='?', help='See final fit of chunk', default=[], const=':100', type=arg2range)
    argopt('-lookguess', nargs='?', help='Show inital model', default=[], const=':100', type=arg2range)
    argopt('-lookpar', nargs='?', help='See parameter of chunk', default=[], const=':100', type=arg2range)
    argopt('-nset', help='index for spectrum', default=':', type=arg2slice)
    argopt('-nexcl', help='Pattern ignore', default=[], type=arg2range)
    argopt('-oset', help='index for order', default=oset, type=arg2slice)
    argopt('-tag', help='Output tag for filename', default='tmp', type=str)
    argopt('-targ', help='Target name requested in simbad for coordinates, proper motion, parallax and absolute RV.', dest='targname')
    argopt('-vg', help='RV guess', default=1., type=float)   # slightly offsetted
    argopt('-?', '-h', '-help', '--help',  help='show this help message and exit', action='help')

    args = parser.parse_args()
    globals().update(vars(args))


def fit_chunk(o, chunk, obsname, targ=None, tpltarg=None):
    ####  observation  ####
    x, w, f, bp, bjd, berv = Spectrum(obsname, o=o, targ=targ)
    i = np.arange(f.size)
    bp[mskatm(w) > 0.1] |= 16

    ####  stellar template  ####
    if tplname:
        w_tpl, f_tpl = Tpl(tplname, o=o, targ=tpltarg)
    else:
        # no template given; model pure iodine
        w_tpl, f_tpl = w, 0*f+1

    lmin = max(w[iset][0], w_tpl[0], w_I2[0])
    lmax = min(w[iset][-1], w_tpl[-1], w_I2[-1])

    # trim the observation to a range valid for the model
    vcut = 100   # [km/s]
    bp[np.log(w) < np.log(lmin)+vcut/c] |= 1
    bp[np.log(w) > np.log(lmax)-vcut/c] |= 1

    if chunks > 1:
        # divide dataset into chunks
        i0 = np.where(bp&1==0)[0][0]   # the first pixel that is not trimmed
        len_ch = int(np.sum(bp&1==0) / chunks)

        bp[i0:i0+chunk*len_ch] |= 32
        bp[i0+(chunk+1)*len_ch:] |= 32

    i_ok = np.where(bp==0)[0]
    x_ok = x[i_ok]
    w_ok = w[i_ok]
    f_ok = f[i_ok]


    modset['icen'] = icen = np.mean(x_ok) + 18   # slight offset, then it converges for CES+TauCet

    # display
    # pre-look raw input
    s = slice(*np.searchsorted(w_I2, [lmin, lmax]))
    s_s = slice(*np.searchsorted(w_tpl, [lmin, lmax]))
    sj = slice(*np.searchsorted(uj_full, np.log([lmin, lmax])))

    # prepare input; convert discrete data to model

    # using the supersampled log(wavelength) space with knot index j
    uj = uj_full[sj]
    iod_j = iod_j_full[sj]

    if demo & 1:
        # plot data, template, and iodine with some scaling
        gplot(w_I2[s], f_I2[s]/np.median(f_I2[s]), 'w l lc 9,', w_tpl[s_s], f_tpl[s_s]/np.median(f_tpl[s_s]), 'w l lc 3,', w, f/np.median(f), 'w lp lc 1 pt 7 ps 0.5')
        pause('demo 1: raw input')


    # convert discrete template into a function
    if tplname:
        S_star = lambda x: np.interp(x, np.log(w_tpl)-berv/c, f_tpl)  # Apply barycentric motion
    else:
        S_star = lambda x: 0*x+1

    IP = IPs[ip]

    # setup the model
    S_mod = model(S_star, uj, iod_j, IP, **modset)

    if demo & 2:
        # plot the IP
        gplot(S_mod.vk, S_mod.IP(S_mod.vk))
        pause('demo 2: default IP')

    if demo & 4:
       # plot again, now the stellar template can be interpolated
       gplot(np.exp(uj), iod_j, S_star(uj)/np.median(S_star(uj)), 'w l lc 9, "" us 1:3 w l lc 3')
       pause('demo 4: stellar template evaluate at uj')


    # an initial parameter set
    v = vg   # a good guess for the stellar RV is needed
    a0 = np.mean(f_ok) / np.mean(S_star(np.log(w_ok))) / np.mean(iod_j)
    a = ag = [a0] + [0]*dega
    b = bg = np.polyfit(x_ok-icen, w_ok, degb)[::-1]
    cg = [0] #* degc
    cc, c0 = (cg, []) if degc else ([], cg)
    s = sg = [Inst.pg['s']]
    if demo:
        a = ag = [a0*1.3] + [0]*dega
        b = bg = [w[0], (w[-1]-w[0])/w.size] # [6128.8833940969, 0.05453566108124]
        b = bg = [*np.polyfit(i[[400,-300]]-icen-10, w[[400,-300]], 1)[::-1]] + [0]*(degb-1)
        s = sg = [s[0]*1.5]

    if ip in ('sg', 'mg'):
        s += [2.]   # exponent of super gaussian 
    elif ip in ('ag',):
        s += [1.]   # skewness parameter (offset to get iterationss)

    if demo & 8:
        # a simple call to the forward model
        # Si_mod = S_mod(x_ok, v=0, a=a, b=b, s=s)
        # show the start guess
        S_mod.show([v,a,b,s], x_ok, f_ok, res=False, dx=0.1)
        pause('demo 8: Smod simple call')

    if demo & 16:
        # A wrapper to fit the continuum
        p_a, _ = S_mod.fit(x_ok, f_ok, a=[a0], v0=vg, b0=bg, s0=sg, res=False, dx=0.1)
        ag[0] = p_a[1][0]
        pause('demo 16: S_a')

    if demo & 32:
        # A wrapper to fit the wavelength solution
        p_b, _ = S_mod.fit(x_ok, f_ok, b=bg[:-1]+[1e-15], v0=vg, a0=a, s0=sg,  res=False, dx=0.1)
        b = p_b[2]
        pause('demo 32: S_b')

    if demo & 64:
        # fit v, a0 and b simultaneously
        p, _ = S_mod.fit(x_ok, f_ok, a=a, b=b, v0=vg, s0=sg, res=False, dx=0.1)
        v, a, b, s = p
        pause('demo 64: S_vab')


    if ip in ('sg', 'ag', 'bnd'):
        # prefit with Gaussian IP
        S_modg = model(S_star, uj, iod_j, IPs['g'], **modset)
        if tplname:
            p, _ = S_modg.fit(x_ok, f_ok, v=vg, a=a, b=b, s=sg[0:1], c=cc, c0=c0)
        else:
            # do not fit for RV
            p, _ = S_modg.fit(x_ok, f_ok, a=a, b=b, s=sg[0:1], c=cc, v0=vg, c0=c0)
        v, a, b, *_ = p
        s = [p[-2][0], *sg[1:]]

    if o in lookguess:
        if demo:
            bg = b
            a = [a0]
        pg = [v, a, bg, s, cc+c0]
        prms = S_mod.show(pg, x_ok, f_ok, res=True, dx=0.1)
        pause('lookguess')


    if kapsig:
        # kappa sigma clipping of outliers
        # len1 = len(f_ok)

        pg = [vg, a, bg, s, cc+c0]
        smod = S_mod(x, *pg)
        resid = (f - smod)
        resid[bp != 0] = 0

        bp[abs(resid) >= (kapsig*np.std(resid))] |= 64
        i_ok = np.where(bp == 0)[0]
        x_ok = x[i_ok]
        w_ok = w[i_ok]
        f_ok = f[i_ok]
        # print("Nr of clipped data   :",len1-len(f_ok),"/",len1)

    if IP == 'bnd':
        # Non parametric fit with band matrix
        # We step through velocity in 100 m/s step. At each step there is linear least square
        # fit for the 2D IP using band matrix.
        S_mod = model_bnd(S_star, uj, iod_j, p[2], **modset)
        opt = {'x': x_ok, 'sig_k': s[0]/1.5/c}
        rr = S_mod.fit(f_ok, 0.1, **opt)
        fx = S_mod(x_ok, 0.1, rr[0])
        ipxj = S_mod.IPxj(rr[0])
        if demo & 2:
            gplot(ipxj, 'matrix w image')

        e_v = np.nan
        if tplname:
            vv = np.arange(-1, 1, 0.1)
            RR = []
            aa = []
            for v in vv:
                rr = S_mod.fit(f_ok, v, **opt)
                RR.append(*rr[1])
                aa.append(rr[0])
                if 1:
                    print(v, rr[1])

            v, e_v, a = SSRstat(vv, RR, plot=1, N=f_ok.size)

        best = S_mod.fit(f_ok, v, **opt)
        fx = S_mod(x_ok, v, best[0])
        pause()
        S_mod.show([v, best[0]], x_ok, f_ok, x2=x_ok)
        res = f_ok - fx
        np.savetxt('res.dat', list(zip(x_ok, res)), fmt="%s")
        prms = np.std(res) / fx.mean() * 100
        if o in look:
            pause()
        return v*1000, e_v*1000, bjd.jd, berv, best[0], np.diag(np.nan*best[0]), prms

    if tplname:
        p, e_p = S_mod.fit(x_ok, f_ok, v, a, bg, s, c=cc, c0=c0, dx=0.1)
#        S = lambda x, v, *abs: S_mod(x, v, abs[:1+dega], abs[1+dega:1+dega+1+degb], abs[1+dega+1+degb:])
#        p, e_p = curve_fit(S, x_ok, f_ok, p0=[v]+a+[*bg]+s, epsfcn=1e-12)
    else:
        # do not fit for velocity
        p, e_p = S_mod.fit(x_ok, f_ok, a=a, b=bg, s=s, c=cc, v0=0, c0=c0, dx=0.1)
        # prepend dummy parameter
#        e_p = np.diag([np.nan, *np.diag(e_p)])

    # overplot flagged and clipped data
    gplot+(x[bp != 0],w[bp != 0], f[bp != 0], 1*(bp[bp != 0] == 64), 'us (lam?$2:$1):3:(int($4)?5:9) w p pt 6 ps 0.5 lc var t "flagged and clipped"')

    # overplot FTS iodine spectrum
    #gplot+(np.exp(uj), iod_j/iod_j.max()*f_ok.max(), 'w l lc 9')

    rvo, e_rvo = 1000*p[0], 1000*np.diag(e_p)[0]**0.5   # convert to m/s
    #prms = S_mod.show([p[0], p[1:1+1+dega], p[2+dega:2+dega+1+degb], p[3+dega+degb:]], x_ok, f_ok, dx=0.1)
    # gplot+(w_tpl[s_s]*(1-berv/c), f_tpl[s_s]*ag, 'w lp lc 4 ps 0.5')
    #gplot+(x_ok, S_star(np.log(np.poly1d(b[::-1])(x_ok))+(v)/c), 'w lp ps 0.5')
    # gplot+(np.exp(S_star.x), S_star.y, 'w lp ps 0.5 lc 7')
    fmod = S_mod(x_ok, *p)
    res = f_ok - fmod
    prms = np.std(res) / np.mean(fmod) * 100
    np.savetxt('res.dat', list(zip(x_ok, res)), fmt="%s")

    if o in look:
        pause('look %s:'% o,  rvo,'+/- %.2f' % e_rvo)  # globals().update(locals())

    sa = tplname is not None
    sb = sa + dega+1
    ss = sb + degb+1
    # error estimation
    # uncertainty in continuum
    xl = np.log(np.poly1d(p[2][::-1])(i-icen))
    Cg = np.poly1d(ag[::-1])(i-icen)      # continuum guess
    Cp = np.poly1d(p[1][::-1])(i-icen)    # best continuum
    X = np.vander(xl, dega+1)[:,::-1].T
    e_Cp = np.einsum('ji,jk,ki->i', X, e_p[sa:sb,sa:sb], X)**0.5
    # uncertainty in wavelength solution
    X = np.vander(xl, degb+1)[:,::-1].T
    lam_g = np.poly1d(bg[::-1])(i-icen)
    lam = np.poly1d(p[2][::-1])(i-icen)
    e_lam = np.einsum('ji,jk,ki->i', X, e_p[sb:ss,sb:ss], X)**0.5
    e_wavesol = np.sum((e_lam/lam*3e8)**-2)**-0.5

    if o in lookpar:
        # compare the wavelength solutions
        #show_model(i, np.poly1d(b[::-1])(i), np.poly1d(bg[::-1])(i), res=True)
        gplot.reset()
        gplot.multiplot("layout 2,2")
        gplot.xlabel('"pixel"').ylabel('"k(x2)"')
        gplot.mxtics().mytics()
        gplot('[][0:]', i, Cg, Cp, e_Cp, 'w l lc 9 t "guess",  "" us 1:3 w l lc 3, "" us 1:($3-$4):($3+$4) w filledcurves fill fs transparent solid 0.2 lc 3 t "1{/Symbol s}" ')
        gplot.xlabel('"pixel"').ylabel('"deviation c * ({/Symbol l} / {/Symbol l}_{guess} - 1) [km/s]"')
        gplot(i, (lam/lam_g-1)*c, ((lam-e_lam)/lam_g-1)*c, ((lam+e_lam)/lam_g-1)*c, 'w l lc 3, "" us 1:3:4 w filledcurves fill fs transparent solid 0.2 lc 3 t "1{/Symbol s}"')
        gplot.xlabel('"[km/s]"').ylabel('"contribution"')
        e_s = e_p[ss,ss]**0.5
        s = p[3]
        gplot(S_mod.vk, S_mod.IP(S_mod.vk, *sg), ' lc 9 ps 0.5 t "IP_{guess}", ', S_mod.vk, S_mod.IP(S_mod.vk,*p[3]), S_mod.IP(S_mod.vk,*[s[0]-e_s, *s[1:]]),  S_mod.IP(S_mod.vk, *[s[0]+e_s, *s[1:]]), 'lc 3 ps 0.5 t "IP", "" us 1:3:4 w filledcurves fill fs transparent solid 0.2 lc 3 t "1{/Symbol s}"')
        gplot.unset('multiplot')
        pause('lookpar')

    return rvo, e_rvo, bjd.jd, berv, p, e_p, prms


obsnames = np.array(sorted(glob.glob(obspath)))[nset]
obsnames = [x for i,x in enumerate(obsnames) if i not in nexcl]
N = len(obsnames)
if not N: pause('no files: ', obspath)

if targname:
    targ = Targ(targname, csv=tag+'.targ.csv').sc

orders = np.r_[oset]
print(orders)

rv = np.nan * np.empty((chunks*len(orders)))
e_rv = np.nan * np.empty((chunks*len(orders)))

rvounit = open(tag+'.rvo.dat', 'w')
parunit = open(tag+'.par.dat', 'w')
# file header
print('BJD', 'RV', 'e_RV', 'BERV', *sum(zip(map("rv{}".format, orders), map("e_rv{}".format, orders)),()), 'filename', file=rvounit)
print('BJD', 'o', *sum(zip(map("p{}".format, range(10)), map("e_p{}".format, range(10))),()), 'prms', file=parunit)

####  FTS  ####

# using the supersampled log(wavelength) space with knot index j

w_I2, f_I2, uj_full, iod_j_full = FTS(ftsname)
mskatm = lambda x: np.interp(x, *np.genfromtxt(viperdir+'lib/mask_vis1.0.dat').T)


T = time.time()
for n,obsname in enumerate(obsnames):
    filename = os.path.basename(obsname)
    print("%2d/%d"% (n+1,N), obsname)
    for i_o, o in enumerate(orders):
        for ch in np.arange(chunks):
            gplot.RV2title = lambda x: gplot.key('title noenhanced "%s (n=%s, o=%s%s)"'% (filename, n+1, o, x))
            gplot.RV2title('')
            rv[i_o*chunks+ch], e_rv[i_o*chunks+ch], bjd,berv, p, e_p, prms = fit_chunk(o, ch, obsname=obsname, targ=targ, tpltarg=targ)
#        try:
#            rv[i_o], e_rv[i_o], bjd,berv, p, e_p  = fit_chunk(o, obsname=obsname)
#        except Exception as e:
#            if repr(e) == 'BdbQuit()':
#               exit()
#            print("Order failed due to:", repr(e))

            print(n+1, o, ch, rv[i_o*chunks+ch], e_rv[i_o*chunks+ch])
            print(bjd, o, ch, *sum(zip(p, np.diag(e_p)),()), prms, file=parunit)
            # store residuals
            os.system('mkdir -p res; touch res.dat')
            os.system('cp res.dat res/%03d_%03d.dat' % (n,o))
            pause()

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

pause('%s done.' % tag)

