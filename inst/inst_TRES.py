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

location = tres = EarthLocation.from_geodetic(lat=31.68094*u.deg, lon=110.8775*u.deg, height=2320*u.m)

oset = '10:40'
#iset = '380:1700'

# convert FHWM resolution to sigma
ip_guess = {'s': 300_000/62_000/ (2*np.sqrt(2*np.log(2))) }   

class observation:

    def __init__(self, filename, targ, *args):
    
        self.filename = filename
        self.hdu = hdu = fits.open(filename, ignore_blank=True)[0]
        self.hdr = hdr = hdu.header

        dateobs = hdr.get('DATE-OBS', 0)
        exptime = hdr.get('EXPTIME', 0)   
        ra = hdr.get('RA', np.nan)                          
        de = hdr.get('DEC', np.nan)
    
        ra = tuple(map(float, ra.split(':')))
        de = tuple(map(float, de.split(':')))
        ra = np.polyval(ra[::-1], 1/60)*15
        de = np.polyval(np.copysign(de[::-1], de[0]), 1/60)

        targdrs = SkyCoord(ra=ra*u.deg, dec=de*u.deg)
    
        if not targ: targ = targdrs
        midtime = Time(dateobs, format='isot', scale='utc') + exptime/2. * u.s
    
        bjd = midtime.tdb
        
        self.bjd, self.targ = bjd, targ

        spec = hdu.data
        gg = readmultispec(filename, reform=True, quiet=True)
        wave = gg['wavelen']
        wave = airtovac(wave)
        
        self.wave_all, self.spec_all = wave_all, spec_all
        
    def Spectrum(self, order):        
        
        if order is not None:
            wave, spec = self.wave_all[order], self.spec_all[order]

        pixel = np.arange(spec.size) 
        err = np.zeros(spec.size)+0.1
        flag_pixel = 1 * np.isnan(spec) # bad pixel map
        #b[spec>1.] |= 4   # large flux, only for normalised spectra, use kapsig instead
        flag_pixel[(5300<wave) & (wave<5343)] |= 256  # only for HARPS s1d template (this order misses)

        return pixel, wave, spec, err, flag_pixel

def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    
    wave, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ) 

    return wave, spec


def FTS(ftsname='None', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)


def write_fits(wtpl_all, tpl_all, e_all, list_files, file_out):

    file_in = list_files[0]

    # copy header from first fits file 
    hdu = fits.open(file_in, ignore_blank=True)[0]
    try:
        # causing problems for some files because of wrong FITS format
        del hdu.header['P.I.']	
    except:
        pass
    
    f = hdu.data

    # write the template data to the file
    for o in range(1, len(f), 1): 
        if o in tpl_all:
            f[o] = tpl_all[o]
        else:
            f[o] = np.ones(len(f[o]))

    hdu.writeto(file_out+'_tpl.model', overwrite=True) 


