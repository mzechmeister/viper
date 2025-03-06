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

location = mcdonald = EarthLocation.from_geodetic(lat=30.671896*u.deg, lon=-104.022143*u.deg, height=2062.46*u.m)
oset = '20:37'

# convert FHWM resolution to sigma
ip_guess = {'s': 300_000/67_000/ (2*np.sqrt(2*np.log(2))) }   

def Spectrum(filename='', order=None, targ=None):
    hdu = fits.open(filename, ignore_blank=True)[0]
    hdr = hdu.header

    try:
        dateobs = hdr['DATE-OBS']+ 'T' + hdr['MIDTIME']
        exptime = 0
    except:
        dateobs = hdr['DATE-OBS']+ 'T' + hdr['UT']
        exptime = hdr.get('EXPTIME')  
        
    ra = hdr.get('RA', np.nan)                          
    de = hdr.get('DEC', np.nan)

    targdrs = SkyCoord(ra=ra, dec=de, unit=(u.hourangle, u.deg))
    if not targ: targ = targdrs
    midtime = Time(dateobs, format='isot', scale='utc') + exptime/2. * u.s
    berv = targ.radial_velocity_correction(obstime=midtime, location=mcdonald)
    berv = berv.to(u.km/u.s).value
    bjd = midtime.tdb

    spec = hdu.data
    gg = readmultispec(filename, reform=True, quiet=True)
    wave = gg['wavelen']
    wave = airtovac(wave)
    if order is not None:
         wave, spec = wave[order], spec[order]

    pixel = np.arange(spec.size) 
    err = np.zeros(spec.size)+0.1
    flag_pixel = 1 * np.isnan(spec) # bad pixel map
    #b[spec>1.] |= 4   # large flux, only for normalised spectra, use kapsig instead

    return pixel, wave, spec, err, flag_pixel, bjd, berv

def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    
    wave, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ)

    return wave, spec


def FTS(ftsname='lib/McDonald/mcd1.fits', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)


def write_fits(wtpl_all, tpl_all, e_all, list_files, file_out):

    file_in = list_files[0]

    # copy header from first fits file 
    hdu = fits.open(file_in, ignore_blank=True)[0]
    f = hdu.data

    # write the template data to the file
    for o in range(1, len(f), 1): 
        if o in tpl_all:
            f[o] = tpl_all[o]
        else:
            f[o] = np.ones(2048)

    hdu.writeto(file_out+'_tpl.model', overwrite=True) 

