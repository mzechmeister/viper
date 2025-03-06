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


oset = '32:82'

# convert FHWM resolution to sigma
ip_guess = {'s': 300_000/67_000/ (2*np.sqrt(2*np.log(2))) }   

def Spectrum(filename='', order=None, targ=None):
    hdu = fits.open(filename, ignore_blank=True)
    hdr = hdu[0].header

    dateobs = hdr.get('DATE-OBS', np.nan)
    ra = hdr.get('RA', np.nan)                        
    de = hdr.get('DEC', np.nan)

    berv = hdr.get('TNG DRS BERV', 0)
    bjd = Time(dateobs, format='isot', scale='utc')

    data = hdu[1].data

    if order is not None:
         wave, spec = data['WAVE'][order-32][::-1]*10, data['FLUX'][order-32][::-1]

    pixel = np.arange(spec.size) 
    err = np.zeros(spec.size)+0.1
    flag_pixel = 1 * np.isnan(spec) # bad pixel map

    return pixel, wave, spec, err, flag_pixel, bjd, berv

def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    
    wave, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ)

    return wave, spec


def FTS(ftsname='', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)


def write_fits(wtpl_all, tpl_all, e_all, list_files, file_out):

    file_in = list_files[0]

    # copy header from first fits file 
    hdu = fits.open(file_in, ignore_blank=True)
    hdr = hdu[0].header
    f = hdu[1].data
    hdr.set('VIPER CAL', 'tell corr')

    # write the template data to the file
    for o in range(32,82,1): 
        if o in tpl_all:
            f['FLUX'][o-32] = tpl_all[o][::-1]
        else:
            f['FLUX'][o-32] = np.ones(2048)

    hdu.writeto(file_out+'_tpl.fits', overwrite=True) 


