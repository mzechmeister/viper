#! /usr/bin/env python3
# Licensed under a GPLv3 style license - see LICENSE

import numpy as np
import os.path
import sys
import os
from datetime import datetime
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

# EXES spectograph at the SOFIA spacecraft
# no location given; height ~13 km
location = exes = EarthLocation.from_geodetic(
    lat=0 * u.deg, lon=0 * u.deg, height=0 * u.m
)

oset = '1:60'

ip_guess = {'s': 300_000/150_000/ (2*np.sqrt(2*np.log(2))) } 

class observation:

    def __init__(self, filename, targ, *args):

        self.filename = filename
        
        if filename.endswith('.fits'):      
            try:
                self.hdu = hdu = fits.open(filename, ignore_blank=True)[0]
             #   self.hdr = hdr = hdu[0].header
                self.data = data = hdu.data
                
                self.format = 1
            
            except:
                self.hdu = hdu = fits.open(filename, ignore_blank=True)
                self.hdr = hdr = hdu[0].header
  
                self.format = 2   
       
        else:

            data = np.genfromtxt(filename, dtype=None, names=True, deletechars='', encoding=None).view(np.recarray)
    
            if 1:
                wave = data.wave
                spec = data.flux
    
        #  wave = np.array([d[0] for d in data])
        # spec = np.array([d[1] for d in data])
    
            try:
                wave, spec = 1e8 /data.wave[::-1], data.flux[::-1]
            
                if 0:
                    blaze = data.blk[::-1]
                    spec /= blaze
           
            except:
                wave = np.array([d[0] for d in data])
                spec = np.array([d[1] for d in data])
                wave = 1e8/wave[::-1]
                spec = spec[::-1]
  
            
            self.wave_all, self.spec_all = wave, spec
            
            self.format = 3
            
      #  dateobs = hdr['DATE-OBS']
        
      #  midtime = Time(dateobs, format='isot', scale='utc') 
      
        dateobs = '2021-06-11T23:31:47.5734'
        midtime = Time(dateobs, format='isot', scale='utc') 
        bjd = midtime.tdb
        berv = 0
        
        self.bjd, self.berv = bjd, berv
  
      
    def Spectrum(self, order):    
    
        match self.format:
            case 1:
                
                wave = self.data[order][0]
                spec = self.data[order][1]
            
                wave = 1e8/wave[::-1]
                spec = spec[::-1]
            
                ps, pe = 170, 710
               # ps, pe = 175, 715
           #     ps, pe = 150, 690
        
                spec = spec[ps:pe]
                wave = wave[ps:pe]
               
                
            case 2:
                data = self.hdu[order].data  
                
                wave = data['wave']
                spec = data['flux'] 
            
            case 3:    
                if "TEXES" in str(self.filename):
                    olen = 1*256	# 1020# 3048
                else:
                    olen = 1020
    
                spec = self.spec_all[order*olen:(1+order)*olen-1]
                wave = self.wave_all[order*olen:(1+order)*olen-1]
                
                self.format = 3
                
        wave = wave[~np.isnan(spec)]
        spec = spec[~np.isnan(spec)]
        pixel = np.arange(spec.size) 
        
        err = np.zeros(spec.size)+0.1
        flag_pixel = 1 * np.isnan(spec) # bad pixel map       

        return pixel, wave, spec, err, flag_pixel


def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''

    if tplname.endswith('.fits'):
        
        try:
            hdu = fits.open(tplname)[0]
            data = hdu.data
        
            wave = data[order][0]
            spec = data[order][1]
            
            wave = 1e8/wave[::-1]
            spec = spec[::-1]
            
        except:
            hdu = fits.open(tplname)[order]
            data = hdu.data
        
            wave = data['wave']
            spec = data['flux']
        
    else:
  
        data = np.genfromtxt(tplname, dtype=None, names=True, deletechars='', encoding=None).view(np.recarray)
    
        try:
            wave, spec = 1e8 /data.wave[::-1], data.flux[::-1]
        except:
            wave = np.array([d[0] for d in data])
            spec = np.array([d[1] for d in data])
            wave = 1e8/wave[::-1]
            spec = spec[::-1]
    
        spec = np.min(spec) +spec

    return wave, spec


def FTS(ftsname='None', dv=100):

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
