#!/usr/bin/env python3

# ./viper.py data/TLS/betgem/BETA_GEM.fits data/TLS/betgem/pepsib.20150409.000.sxt.awl.all6
# ./viper.py data/TLS/hd189733/TV00001.fits data/TLS/Deconv/HD189733.model

import argparse

import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit
from astropy.io import fits

from gplot import *
gplot.tmp = '$'

from inst.inst_TLS import Spectrum, Tpl, FTS
from model import model, IP, show_model

c = 3e5   # [km/s] speed of light


o = 33; lmin = 6120; lmax = 6250
o = 18; lmin = 5240; lmax = 5390

dirname = r''
ftsname = dirname + 'lib/TLS/FTS/TLS_I2_FTS.fits'
obsname = dirname + 'data/TLS/betgem/BETA_GEM.fits'
tplname = dirname + 'data/TLS/betgem/pepsib.20150409.000.sxt.awl.all6'
obsname = dirname + 'data/TLS/hd189733/TV00001.fits'
tplname = dirname + 'data/TLS/Deconv/HD189733.model'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='VIPER - velocity and IP Estimator', add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    argopt = parser.add_argument   # function short cut
    argopt('obsname', help='Filename of observation', default='data/TLS/betgem/BETA_GEM.fits', type=str)
    argopt('tpl', help='Filename of template', default='data/TLS/betgem/pepsib.20150409.000.sxt.awl.all6', type=str)
    argopt('-o', help='index for order', default=18, type=int)
    
    args = parser.parse_args()
    globals().update(vars(args))

####  FTS  ####
w_I2, f_I2 = FTS()

# display
s = slice(*np.searchsorted(w_I2, [lmin, lmax]))
gplot(w_I2[s], f_I2[s], 'w l lc 9')

####  data TLS  ####
w_i, f_i = Spectrum(obsname, o=o)
i = np.arange(f_i.size)

####  stellar template  ####
w_tpl, f_tpl = Tpl(tplname, o=o)


lmin = max(w_tpl[0], w_i[0], w_I2[0])
lmax = min(w_tpl[-1], w_i[-1], w_I2[-1])

s = slice(*np.searchsorted(w_I2, [lmin, lmax]))
s_s = slice(*np.searchsorted(w_tpl, [lmin, lmax]))
gplot(w_I2[s], f_I2[s], 'w l lc 9,', w_tpl[s_s], f_tpl[s_s], 'w l lc 3')


# pre-look data
gplot(w_I2[s], f_I2[s], 'w l lc 9,', w_tpl[s_s], f_tpl[s_s], 'w l lc 3,', w_i, f_i, 'w lp lc 1 pt 7 ps 0.5')
gplot(w_I2[s], f_I2[s]/1.18, 'w l lc 9,', w_tpl[s_s]*(1+12/c), f_tpl[s_s], 'w l lc 3,', w_i, f_i/1.04, 'w lp lc 1 pt 7 ps 0.5')


# prepare input; convert discrete data to model

# define a supersampled log(wavelength) space with knot index j
xj = np.linspace(np.log(lmin)+100/c, np.log(lmax)-100/c, w_I2[s].size)  # reduce range by 100 km/s
iod_j = interpolate.interp1d(np.log(w_I2), f_I2)(xj)


# convert PESPI data into a function
S_star = interpolate.interp1d(np.log(w_tpl), f_tpl)

S_mod = model(S_star, xj, iod_j, IP)

# plot
gplot(S_mod.vk, S_mod.IP(S_mod.vk))

# plot again, now the stellar template can be interpolated
gplot(np.exp(xj), iod_j, S_star(xj), 'w l lc 9, "" us 1:3 w l lc 3')

#gplot(np.exp(xj), iod_j, S_star(xj+3/c), 'w l lc 9 t "iodine", "" us 1:3 w l lc 3 t "template + 3 km/s",', np.exp(xj_eff), S_eff(xj_eff), 'w l lc 1 t "IP x (tpl*I2)"')

# Now wavelength solution

# mapping between pixel and wavelength

#lam(x) = b0 + b_1 * x + b_2 * x**2

lam = np.poly1d([w_i[0], (w_i[-1]-w_i[0])/w_i.size][::-1])
s_obs = slice(*np.searchsorted(np.log(w_i), [xj[0]+100/c, xj[-1]-100/c]))


# well, we see the wavelength solution can be improved

#gplot(i[s_obs], S_eff(np.log(lam(i[s_obs]))), 'w l,', i, f_i, 'w lp pt 7 ps 0.5 lc 3')

# and let's plot the observation against wavelength

#gplot(np.exp(xj_eff), S_eff(xj_eff), 'w l,', lam(i), f_i, 'w lp pt 7 ps 0.5 lc 3')

v=0
a = [0.96]
b = [w_i[0], (w_i[-1]-w_i[0])/w_i.size] # [6128.8833940969, 0.05453566108124]
s = 2.5


# a simple call to the forward model
Si_mod = S_mod(i[s_obs], v=0, a=[1], b=b, s=s)

#gplot(i, Si_mod, 'w l t "S(i)",', i, f_i, 'w lp pt 7 ps 0.5 lc 3 t "S_i"')
show_model(i[s_obs], f_i[s_obs], Si_mod, res=False)

# A wrapper to fit the continuum
S_a = lambda x, a0: S_mod(x, v, [a0], b, s)

a, e_a = curve_fit(S_a, i[s_obs], f_i[s_obs])

show_model(i[s_obs], f_i[s_obs], S_a(i[s_obs],*a), res=False)

# A wrapper to fit the wavelength solution
S_b = lambda x, b0,b1,b2,b3: S_mod(x, v, a, [b0,b1,b2,b3], s)

v = -2.   # a good guess for the stellar RV is needed
bg = np.polyfit(i[s_obs], w_i[s_obs], 3)[::-1]
b, e_b = curve_fit(S_b, i[s_obs], f_i[s_obs], p0=bg)
bg1 = b*1

#show_model(i[s_obs], f_i[s_obs], S_b(i[s_obs], *bg))
show_model(i[s_obs], f_i[s_obs], S_b(i[s_obs], *b))
gplot+(i[s_obs], S_star(np.log(np.poly1d(b[::-1])(i[s_obs]))+(v)/c), 'w lp ps 0.5')

# compare the wavelength solutions
show_model(i, np.poly1d(b[::-1])(i), np.poly1d(bg[::-1])(i), res=True)

# fit a and b simulatenously

S_vab = lambda x, v, a, b0,b1,b2,b3: S_mod(x, v, [a], [b0,b1,b2,b3], 2.2)
p, e_p = curve_fit(S_vab, i[s_obs], f_i[s_obs], p0=[v, 1, *bg])
show_model(i[s_obs], f_i[s_obs], S_vab(i[s_obs], *p))
p1 = 1*p

S_vabs = lambda x, v, a, b0,b1,b2,b3, s: S_mod(x, v, [a], [b0,b1,b2,b3], s)
p, e_p = curve_fit(S_vabs, i[s_obs], f_i[s_obs], p0=[*p1, 2.2])
show_model(i[s_obs], f_i[s_obs], S_vabs(i[s_obs], *p))


#show_model(i[s_obs], f_i[s_obs], S_b(i[s_obs], *bg))
show_model(i[s_obs], f_i[s_obs], S_vabs(i[s_obs], *p))
gplot+(i[s_obs], S_star(np.log(np.poly1d(b[::-1])(i[s_obs]))+(v)/c), 'w lp ps 0.5')


print('Done.') 
