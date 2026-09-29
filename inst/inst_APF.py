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

path = os.path.join(os.path.dirname(__file__), '..') + "/lib/APF/"

# Automated Planet Finder @Lick, Mount Hamilton California
location = APF = EarthLocation.from_geodetic(lat=37.343333*u.deg, lon=-121.6366667*u.deg, height=1290*u.m)

oset = '32:44'

# convert FHWM resolution to sigma
ip_guess = {'s': 300_000/100_000/ (2*np.sqrt(2*np.log(2))) }   

class observation:

    def __init__(self, filename, targ, *args):
    
        self.filename = filename
        self.hdu = hdu = fits.open(filename, ignore_blank=True)[0]
        self.hdr = hdr = hdu.header

        dateobs = hdr.get('THEMIDPT', None)
        exptime = hdr.get('EXPTIME', 0)   
        ra = hdr.get('RA', np.nan)                         
        de = hdr.get('DEC', np.nan)
     
        ra = tuple(map(float, ra.split(':')))
        de = tuple(map(float, de.split(':')))
        ra = np.polyval(ra[::-1], 1/60)*15
        de = np.polyval(np.copysign(de[::-1], de[0]), 1/60)
 
        targdrs = SkyCoord(ra=ra*u.deg, dec=de*u.deg)
        if not targ: targ = targdrs
        midtime = Time(dateobs, format='isot', scale='utc') #+ exptime/2. * u.s
    
        bjd = midtime.tdb  
        self.bjd, self.targ = bjd, targ
    
        # APF uses same start wavelengths for all observations
        # final wavelength calibration is done via the iodine cell
        hdu_wave = fits.open(path+'apf_wave_bj2.fits', ignore_blank=True)[0]
        wave_all = hdu_wave.data  
        spec_all = hdu.data
        
        self.wave_all, self.spec_all = wave_all, spec_all
        
    def Spectrum(self, order):    
 
        if order is not None:
            wave, spec = self.wave_all[order], self.spec_all[order]
            
        wave = airtovac(wave)

        pixel = np.arange(spec.size) 
        err = np.zeros(spec.size)+0.1
        flag_pixel = 1 * np.isnan(spec) # bad pixel map

        return pixel, wave, spec, err, flag_pixel

def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    
    wave, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ) 
    
    return wave, spec

def FTS(ftsname='lib/APF/nist2apf.fits', dv=100):
    
    return resample(*FTSfits(ftsname), dv=dv)

def write_fits(wtpl_all, tpl_all, e_all, list_files, file_out):

    file_in = list_files[0]

    # copy header from first fits file
    hdu = fits.open(file_in, ignore_blank=True)
    hdr = hdu[0].header

    hdr = fits.Header()

    hdu0 = fits.PrimaryHDU()
    
    # write data to the file
    table_all = [hdu0]
    
    o_max = np.max(list(tpl_all.keys()))
    
    for order in np.arange(0, o_max+1):
 #   for order in tpl_all:
        if order in list(tpl_all.keys()):
            c1 = fits.Column(name='wave', array=wtpl_all[order], format='F')
            c2 = fits.Column(name='flux', array=tpl_all[order], format='F')
        else:
            c1 = fits.Column(name='wave', array=wtpl_all[o_max]*0, format='F')
            c2 = fits.Column(name='flux', array=tpl_all[o_max]*0, format='F')

        table_all.append(fits.BinTableHDU.from_columns([c1, c2]))

    hdul = fits.HDUList(table_all)
    
    hdul.writeto(file_out+'_tpl.fits', overwrite=True)
    hdul.close()


