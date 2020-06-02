#!/usr/bin/env python3

import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit
from astropy.io import fits

from gplot import *
gplot.tmp = '$'


o = 33; lmin = 6120; lmax = 6250
#o = 18; lmin = 6120; lmax = 6250


root_FTS = r'lib/TLS/'

c = 3e5   # [km/s] speed of light

#########################  FTS  ########################################

hdu_I2 = fits.open(root_FTS+'/FTS/TLS_I2_FTS.fits')[0]

f_I2 = hdu_I2.data[::-1]
h = hdu_I2.header
w_I2 = h['CRVAL1'] + h['CDELT1'] * (np.arange(f_I2.size) + 1. - h['CRPIX1'])
w_I2 = 1e8 / w_I2[::-1]   # convert wavenumbers to wavelength [angstrom]

# display
s = slice(*np.searchsorted(w_I2, [lmin, lmax]))
gplot(w_I2[s], f_I2[s], 'w l lc 9')


#####  stellar template   ####

hdu = fits.open('data/TLS/other/pepsib.20150409.000.sxt.awl.all6')

w_tpl = hdu[1].data.field('Arg')
f_tpl = hdu[1].data.field('Fun')

s_s = slice(*np.searchsorted(w_tpl, [lmin, lmax]))
gplot(w_I2[s], f_I2[s], 'w l lc 9,', w_tpl[s_s], f_tpl[s_s], 'w l lc 3')



#### data TLS
hdu = fits.open('data/TLS/other/pepsib.20150409.000.sxt.awl.all6')

from inst.inst_TLS import Spectrum
w, f = Spectrum('data/TLS/other/BETA_GEM.fits')
#hdu = fits.open(root_FTS+'BETA_GEM.fits')[0]
#hdu = fits.open('lib/TLS/observations/SPECTRA/TV00007.fits')[0]
#f_i = hdu.data[33]
w_i = w[33]
f_i = f[33]
i = np.arange(f_i.size)
#w_i = 6128.8833940969 + 0.05453566108124*np.arange(f_i.size)  # guess


# pre-look data
gplot(w_I2[s], f_I2[s], 'w l lc 9,', w_tpl[s_s], f_tpl[s_s], 'w l lc 3,', w_i, f_i, 'w lp lc 1 pt 7 ps 0.5')
gplot(w_I2[s], f_I2[s]/1.18, 'w l lc 9,', w_tpl[s_s]*(1+12/c), f_tpl[s_s], 'w l lc 3,', w_i, f_i/1.04, 'w lp lc 1 pt 7 ps 0.5')


# prepare input; convert discrete data to model

# define a supersampled log(wavelength) space with knot index j
xj = np.linspace(np.log(lmin), np.log(lmax), w_I2[s].size)
iod_j = interpolate.interp1d(np.log(w_I2), f_I2/1.18)(xj)

dx = xj[1] - xj[0]  # sampling in uniform resampled Iod
print("sampling [km/s]:", dx*c)

# convert PESPI data into a function
S_star = interpolate.interp1d(np.log(w_tpl), f_tpl)


# IP sampling in velocity space
# index k for IP space
vk = np.arange(-50,50+1) * dx * c
IP_k = np.exp(-(vk/1.5)**2)  # Gauss IP
IP_k /= IP_k.sum()           # normalise IP

# plot
gplot(vk, IP_k)

# plot again, now the stellar template can be interpolated
gplot(np.exp(xj), iod_j, S_star(xj), 'w l lc 9, "" us 1:3 w l lc 3')

# convolving with IP will reduce the valid wavelength range
cut = int(IP_k.size/2)
xj_eff = xj[cut:-cut]

# forward model

def _S_eff(v):
    #S_eff(ln(lam)) = IP(v) x ((a0+a1*ln(lam))*S_PEPSI(ln(lam)+v_star/c) * G_I2(ln(lam)))
    #S_eff(ln(lam)) = IP(v) x ((A+ln(lan) )*S_PEPSI(ln(lam)+v_star/c) * G_I2(ln(lam))) 
    # discrete supersampled effective spectrum
    Sj_eff = np.convolve(IP_k, S_star(xj+v/c) * iod_j, mode='valid')
    
    # return continous supersampled effective spectrum
    return interpolate.interp1d(xj_eff, Sj_eff)


# a forward model for RV shift 3 km/s
S_eff = _S_eff(v=3)


gplot(np.exp(xj), iod_j, S_star(xj+3/c), 'w l lc 9 t "iodine", "" us 1:3 w l lc 3 t "template + 3 km/s (PEPSI)",', np.exp(xj_eff), S_eff(xj_eff), 'w l lc 1 t "IP x (tpl*I2)"')

# Now wavelength solution

# mapping between pixel and wavelength

#lam(x) = b0 + b_1 * x + b_2 * x**2
lam = np.poly1d([6128.8833940969, 0.05453566108124][::-1])

# well, we see the wavelength solution can be improved

gplot(i, S_eff(np.log(lam(i))), 'w l,', i, f_i, 'w lp pt 7 ps 0.5 lc 3')

# and let's plot the observation against wavelength

gplot(np.exp(xj_eff), S_eff(xj_eff), 'w l,', lam(i), f_i, 'w lp pt 7 ps 0.5 lc 3')

v=0
a = [0.96]
b = [6128.8833940969, 0.05453566108124]

def S_mod(i, v, a, b, IP_k):
    '''A forward model'''
    # wavelength solution
    xi = np.log(np.poly1d(b[::-1])(i))
    # IP convolution
    Sj_eff = np.convolve(IP_k, S_star(xj+v/c) * iod_j, mode='valid')
    # sampling to pixel
    Si_eff = interpolate.interp1d(xj_eff, Sj_eff)(xi)
    # flux normalisation
    Si_mod = np.poly1d(a[::-1])(xi) * Si_eff
    return Si_mod


def show_model(x, y, ymod, res=True):
    gplot(x, y, ymod, 'w lp pt 7 ps 0.5 t "S_i",',
          '"" us 1:3w lp pt 6 ps 0.5 lc 3 t "S(i)"')
    if res:
        gplot.mxtics().mytics().my2tics()
        # overplot residuals
        gplot.y2range('[-0.2:2]').ytics('nomirr').y2tics()
        gplot+(x, y-ymod, "w p pt 7 ps 0.5 lc 1 axis x1y2 t 'res', 0 lc 3axis x1y2")

# a simple call to the forward model
Si_mod = S_mod(i, v=0, a=[1], b=b, IP_k=IP_k)

gplot(i, Si_mod, 'w l t "S(i)",', i, f_i, 'w lp pt 7 ps 0.5 lc 3 t "S_i"')
show_model(i, f_i, Si_mod, res=False)

# A wrapper to fit the continuum
S_a = lambda x, a0: S_mod(x, v, [a0], b, IP_k)

a, e_a = curve_fit(S_a, i, f_i)

show_model(i, f_i, S_a(i,*a), res=False)

# A wrapper to fit the wavelength solution
S_b = lambda x, b0,b1,b2,b3: S_mod(x, v, a, [b0,b1,b2,b3], IP_k)

v = -12   # a good guess for the stellar RV is needed
bg = np.polyfit(i, w_i, 3)[::-1]
b, e_b = curve_fit(S_b, i, f_i, p0=bg)

show_model(i, f_i, S_b(i,*bg))

# compare the wavelength solutions
show_model(i, np.poly1d(b[::-1])(i), np.poly1d(bg[::-1])(i), res=True)



 
