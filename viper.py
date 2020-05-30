#!/usr/bin/env python3

import numpy as np
from scipy import interpolate
from astropy.io import fits

from gplot import *
gplot.tmp = '$'

lmin = 6120
lmax = 6250

root_FTS = r'lib/TLS/other/'


#########################  FTS  ########################################

hdu_I2 = fits.open(root_FTS+'/../FTS/TLS_I2_FTS.fits')[0]

f_I2 = hdu_I2.data[::-1]
h = hdu_I2.header
w_I2 = h['CRVAL1'] + h['CDELT1'] * np.arange(f_I2.size)   # re-check conversion
w_I2 = 1e8 / w_I2[::-1]   # scale convertion from wavenumber to wavelength (angstrom)

# display
s = slice(*np.searchsorted(w_I2, [lmin, lmax]))
gplot(w_I2[s], f_I2[s], 'w l lc 9')


#####  stellar template   ####

hdu = fits.open(root_FTS+'pepsib.20150409.000.sxt.awl.all6')

w_tpl = hdu[1].data.field('Arg')
f_tpl = hdu[1].data.field('Fun')

s_s = slice(*np.searchsorted(w_tpl, [lmin, lmax]))
gplot(w_I2[s], f_I2[s], 'w l lc 9,', w_tpl[s_s], f_tpl[s_s], 'w l lc 3')



#### data TLS
#34 92 2 6128.8833940969 0.05453566108124 2048 0.'
#WAT2_102= ' 1165.31 1187.31 1. 0. 2 6 1. 2048. 6185.83841636954 55.972580248164'

hdu = fits.open(root_FTS+'BETA_GEM.fits')[0]
#hdu = fits.open('lib/TLS/observations/SPECTRA/TV00007.fits')[0]
f_i = hdu.data[33]
i = np.arange(f_i.size)
w_i = 6128.8833940969 + 0.05453566108124*np.arange(f_i.size)  # guess


# pre-look data
gplot(w_I2[s], f_I2[s], 'w l lc 9,', w_tpl[s_s], f_tpl[s_s], 'w l lc 3,',w_i,f_i, 'w lp lc 1 pt 7 ps 0.5')
gplot(w_I2[s], f_I2[s]/1.18, 'w l lc 9,', w_tpl[s_s]*(1+12/3e5), f_tpl[s_s], 'w l lc 3,',w_i,f_i/1.04, 'w lp lc 1 pt 7 ps 0.5')


# prepare input; convert discrete data to model

# define a supersampled log(wavelength) space with knot index k
xk = np.linspace(np.log(lmin), np.log(lmax), w_I2[s].size)
iod_k = interpolate.interp1d(np.log(w_I2), f_I2/1.18)(xk)

dx = xk[1] - xk[0]  # sampling in uniform resampled Iod
print("sampling [km/s]:", dx*3e5)

# convert PESPI data into a function
S_star = interpolate.interp1d(np.log(w_tpl), f_tpl)


# IP sampling in velocity space
vl = np.arange(-50,50+1) * dx * 3e5
IP_l = np.exp(-(vl/1.5)**2)  # Gauss IP
IP_l /= IP_l.sum()           # normalise IP

# plot
gplot(vl, IP_l)

# plot again, now the stellar template can be interpolated
gplot(np.exp(xk), iod_k, S_star(xk), 'w l lc 9, "" us 1:3 w l lc 3')

# convolving with IP will reduce the valid wavelength range
cut = int(IP_l.size/2)
xk_eff = xk[cut:-cut]

# forward model

def _S_eff(v):
    #S_eff(ln(lam)) = IP(v) x ((a0+a1*ln(lam))*S_PEPSI(ln(lam)+v_star/c) * G_I2(ln(lam)))
    #S_eff(ln(lam)) = IP(v) x ((A+ln(lan) )*S_PEPSI(ln(lam)+v_star/c) * G_I2(ln(lam))) 
    # discrete supersampled effective spectrum
    Sk_eff = np.convolve(IP_l, S_star(xk+v/3e5) * iod_k, mode='valid')
    
    # return continous supersampled effective spectrum
    return interpolate.interp1d(xk_eff, Sk_eff)


# a forward model for RV shift 3 km/s
S_eff = _S_eff(v=3)


gplot(np.exp(xk), iod_k, S_star(xk+3/3e5), 'w l lc 9 t "iodine I", "" us 1:3 w l lc 3 t "template + 3 km/s (PEPSI)",', np.exp(xk_eff), S_eff(xk_eff), 'w l lc 1 t "IP x (tpl*I2)"' )

# Now wavelength solution

# mapping between pixel and wavelength

#lam(x) = b0 + b_1 * x + b_2 * x**2
lam = np.poly1d([6128.8833940969, 0.05453566108124][::-1])

# well, we see the wavelength solution can be improved

gplot(i, S_eff(np.log(lam(i))), 'w l,', i, f_i, 'w lp pt 7 ps 0.5 lc 3')

# ok lets plot the observation against wavelength,

gplot(np.exp(xk_eff), S_eff(xk_eff), 'w l,', lam(i), f_i, 'w lp pt 7 ps 0.5 lc 3')

