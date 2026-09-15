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

oset = '1:100'
iset = '10:-10'

# convert FHWM resolution to sigma
ip_guess = {'s': 300_000/80_000/ (2*np.sqrt(2*np.log(2))) }   

class observation:

    def __init__(self, filename, targ, *args):
    
        self.filename = filename
        self.hdu = hdu = fits.open(filename, ignore_blank=True)
        self.hdr = hdr = hdu[0].header

        dateobs = hdr.get('DATE-OBS', np.nan)
        exptime = hdr.get('EXP', 0)   
        ra = hdr.get('RA', np.nan)                        
        de = hdr.get('DEC', np.nan)
        tel_lat = hdr.get('TEL GEOLAT', np.nan)
        tel_lon = hdr.get('TEL GEOLON', np.nan)
        tel_hei = hdr.get('TEL GEOELEV', np.nan)

        if len(dateobs) < 12:
            print('WARNING: Incorrect time information in FITS header. This will lead to a wrong barycentric correction.')

        targdrs = SkyCoord(ra=ra*u.hour, dec=de*u.deg)
        if not targ: targ = targdrs
        midtime = Time(dateobs, format='isot', scale='utc') + exptime/2. * u.s
    
        self.location = location = EarthLocation.from_geodetic(lat=tel_lat*u.deg, lon=tel_lon*u.deg, height=tel_hei*u.m)
    
        berv = targ.radial_velocity_correction(obstime=midtime, location=location)
        berv = berv.to(u.km/u.s).value
        bjd = midtime.tdb
    
        self.bjd, self.targ, self.berv = bjd, targ, berv
    
    def Spectrum(self, order):
        if order is not None:
            wave = self.hdu[order].data['wave']
            spec = self.hdu[order].data['flux']
         
            try:
                err = hdu[order].data['error']
            except:
                err = spec*0+0.1
         
   # wave = airtovac(wave)

        pixel = np.arange(spec.size) 
        flag_pixel = 1 * np.isnan(spec) # bad pixel map

        return pixel, wave, spec, err, flag_pixel
    
def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    wave, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ) 

    return wave, spec


def FTS(ftsname='', dv=100):

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
            f[o] = np.ones(len(f[o]))

    hdu.writeto(file_out+'_tpl.model', overwrite=True) 


