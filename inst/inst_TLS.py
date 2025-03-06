#! /usr/bin/env python3
# Licensed under a GPLv3 style license - see LICENSE

import numpy as np
import sys
import os
from astropy.io import fits
from astropy.time import Time
import datetime
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.constants import c

from .template import read_tpl
from .readmultispec import readmultispec
from .airtovac import airtovac

from .FTS_resample import resample, FTSfits

# see https://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py

location = tls = EarthLocation.from_geodetic(lat=50.980111*u.deg, lon=11.711167*u.deg, height=342*u.m)

oset = '18:30'
iset = '380:1700'

# convert FHWM resolution to sigma
ip_guess = {'s': 300_000/67_000/ (2*np.sqrt(2*np.log(2))) }   

def Spectrum(filename='data/TLS/other/BETA_GEM.fits', order=None, targ=None):
    hdu = fits.open(filename, ignore_blank=True)[0]
    hdr = hdu.header

    dateobs = hdr.get('DATE-OBS', hdr.get('FRAME'))
    exptime = hdr.get('EXP_TIME', hdr.get('EXPOSURE'))   # 20211018_guenther_TCEcell_0063.fits EXPOSURE (no exptime)
    ra = hdr.get('RA', np.nan)                          # 20211018_guenther_TCEcell_0184.fits no RA
    de = hdr.get('DEC', np.nan)
    
    if len(dateobs) < 12:
        # in old data, date and time of the observation are stored separately
        try: 
            tm = hdr.get('TM-START')
            dateobs = str(dateobs) + 'T' + str(datetime.timedelta(seconds=tm))
        except:
            print('WARNING: Incorrect time information in FITS header. This will lead to a wrong barycentric correction.')

    # 2021 the ANCHOR CCD was installed at TLS, which writes the end-time of 
    # the observation instead of the start-time
    # -> therefore the exptime has to be subracted instead of added  
    if dateobs > '2021-03-15-00T00:00:00':
        exptime *= -1.

    targdrs = SkyCoord(ra=ra*u.hour, dec=de*u.deg)
    if not targ: targ = targdrs
    midtime = Time(dateobs, format='isot', scale='utc') + exptime/2. * u.s
    
    berv = targ.radial_velocity_correction(obstime=midtime, location=tls)
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
    flag_pixel[(5300<wave) & (wave<5343)] |= 256  # only for HARPS s1d template (this order misses)
    # TLS spectra have a kink in continuum  at about 1700
    # Also the deconv could have a bad wavelength solution.

    return pixel, wave, spec, err, flag_pixel, bjd, berv

def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    
    wave, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ) 

    return wave, spec


def FTS(ftsname='lib/TLS/FTS/TLS_I2_FTS.fits', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)


def write_fits(wtpl_all, tpl_all, e_all, list_files, file_out):

    file_in = list_files[0]

    # copy header from first fits file 
    hdu = fits.open(file_in, ignore_blank=True)[0]
    try:
        # causing problems for some files because of wrong FITS format
        del hdu.header['UT']	
    except:
        pass
    
    f = hdu.data

    # write the template data to the file
    for o in range(1, len(f), 1): 
        if o in tpl_all:
            f[o] = tpl_all[o]
        else:
            f[o] = np.ones(2048)

    hdu.writeto(file_out+'_tpl.model', overwrite=True) 


