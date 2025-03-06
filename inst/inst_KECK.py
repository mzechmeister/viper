#! /usr/bin/env python3
# Licensed under a GPLv3 style license - see LICENSE

import numpy as np
import sys
import os
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.constants import c

from .template import read_tpl
from .readmultispec import readmultispec
from .airtovac import airtovac

from .FTS_resample import resample, FTSfits

# see https://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py

location = keck = EarthLocation.of_site('Keck Observatory')

oset = '5:16'
iset = '500:2000'

ip_guess = {'s': 300_000/87_000/ (2*np.sqrt(2*np.log(2))) }   # convert FHWM resolution to sigma

def Spectrum(filename, order=None, targ=None):
    if order is not None:
         filename = filename.replace('_flux.fits.gz', '')[:-2]+"%02i_flux.fits.gz" % order

    hdu = fits.open(filename, ignore_blank=True)
    hdr = hdu[0].header

    dateobs = hdr['DATE_BEG']   # ~ DATE-OBS+UTC
    exptime = hdr['EXPTIME']
    ra = hdr['RA']
    de = hdr['DEC']
    iod_in = hdr['IODIN']

    targdrs = SkyCoord(ra=ra, dec=de, unit=(u.deg,u.deg))
    if not targ: targ = targdrs
    midtime = Time(dateobs, format='isot', scale='utc') + exptime * u.s
    berv = targ.radial_velocity_correction(obstime=midtime, location=keck)
    berv = berv.to(u.km/u.s).value
    bjd = midtime.tdb

    d = hdu[1].data
    w, f, e = d['wave'], d['Flux'], d['Error']
    w *= (1-(0+berv)/3e5)
    #from pause import pause;   # pause()

    x = np.arange(f.size) 
    #f = np.array([*fits.open('/home/astro115/carmenes/data/KECK/red/HI.20110118.42499_x.fits')[0].data]+[0]*56)
    #f = np.array([*fits.open('/home/astro115/carmenes/data/KECK/red/HI.20060412.39465_x.fits')[0].data]+[0]*56)
    #pause(); from gplot import gplot
    #f = np.array([*fits.open('/home/astro115/carmenes/data/KECK/red/hd189733_test//HI.20060821.19898_x2.fits')[1].data[o-1]]+[0]*56)
    #f[[1351,1352,1662,1663]] = np.nan
    b = 1 * np.isnan(f) # bad pixel map
    #b[f>1.5] |= 2 # large flux
    b[(5300<w) & (w<5343)] |= 256  # only for HARPS s1d template (this order misses)
    # TLS spectra have a kink in continuum  at about 1700
    # Also the deconv could have a bad wavelength solution.


    return x, w, f, e, b, bjd, berv
    '''
     hdr = hdu.header

    dateobs = '2020-08-11'
    exptime = hdr['EXP_TIME']
    ra = hdr['RA']
    de = hdr['DEC']

    targdrs = SkyCoord(ra=ra*u.hour, dec=de*u.deg)
    if not targ: targ = targdrs
    midtime = Time(dateobs, format='isot', scale='utc') + exptime * u.s
    berv = targ.radial_velocity_correction(obstime=midtime, location=tls)
    berv = berv.to(u.km/u.s).value
    bjd = midtime.tdb
    '''

    if filename.endswith('.dat'):
        bjd = 0.
        berv = 0.
        w, f = np.loadtxt(filename).T

        if 0 and order is not None:
            w, f = w[order], f[order]

        x = np.arange(f.size) 
        e = np.ones(f.size)
        b = 1 * np.isnan(f) # bad pixel map
        b[f>1e6] |= 2 # large flux
        b[(5300>w) | (w>5343)] |= 4  # only for HARPS s1d template (this order misses)

    return x, w, f, e, b, bjd, berv

def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    
    wave, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ)

    return w, f


def FTS(ftsname='lib/KECK/fts_keck50_not_normalized.dat', dv=100):
    w, f = np.loadtxt(ftsname).T
    return resample(w, f, dv=dv)
'''
def FTS(ftsname='lib/TLS/FTS/TLS_I2_FTS.fits', dv=100):
    w, f = FTSfits(ftsname)
    return resample(w*(1-0/3e5),f, dv=dv)
'''


