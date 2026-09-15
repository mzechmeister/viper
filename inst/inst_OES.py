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

location = oes = EarthLocation.from_geodetic(lat=49.91056*u.deg, lon=14.78361*u.deg, height=528*u.m)

oset = '1:30'

ip_guess = {'s': 300_000/67_000/ (2*np.sqrt(2*np.log(2))) }   # convert FHWM resolution to sigma

class observation:

    def __init__(self, filename, targ, *args):
    
        self.filename = filename
        self.hdu = hdu = fits.open(filename, ignore_blank=True)[0]
        self.hdr = hdr = hdu.header

        dateobs = hdr['DATE-OBS']+ 'T' + hdr['UT']
        exptime = hdr['EXPTIME']
    
        ra = hdr.get('RA', np.nan)                          
        de = hdr.get('DEC', np.nan)
    
        ra = tuple(map(float, ra.split(':')))
        de = tuple(map(float, de.split(':')))
        ra = np.polyval(ra[::-1], 1/60)*15
        de = np.polyval(np.copysign(de[::-1], de[0]), 1/60)

        targdrs = SkyCoord(ra=ra*u.deg, dec=de*u.deg)
        if not targ: targ = targdrs

        midtime = Time(dateobs, format='isot', scale='utc') + exptime * u.s
        bjd = midtime.tdb
        
        self.bjd, self.targ = bjd, targ
        
        # read in the complete data set
        spec_all = self.hdu.data
        spec_all /= np.nanmean(spec_all)
        gg = readmultispec(self.filename, reform=True, quiet=True)
        wave_all = gg['wavelen']
        wave_all = airtovac(wave_all)
        
        self.wave_all, self.spec_all = wave_all, spec_all
        
    def Spectrum(self, order):
        
        if order is not None:
            wave, spec = self.wave_all[order], self.spec_all[order]

        pixel = np.arange(spec.size) 
        err = np.ones(spec.size)*0.1
        flag_pixel = 1 * np.isnan(spec) # bad pixel map

        return pixel, wave, spec, err, flag_pixel


def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    
    wave, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ)
    
    return wave, spec


def FTS(ftsname='lib/oes.fits', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)


def write_fits(wtpl_all, tpl_all, e_all, list_files, file_out):

    file_in = list_files[0]

    # copy header from first fits file 
    hdu = fits.open(file_in, ignore_blank=True)[0]
    f = hdu.data

    # write the template data to the file
    for order in range(1, len(f), 1): 
        if order in tpl_all:
            f[order] = tpl_all[order]
        else:
            f[order] = np.ones(len(f[order]))

    hdu.writeto(file_out+'_tpl.model', overwrite=True)  
  #  hdu.close()  

