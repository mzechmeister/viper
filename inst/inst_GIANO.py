#! /usr/bin/env python3

## viper - Velocity and IP Estimator
## Copyright (C) Mathias Zechmeister and Jana Koehler
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along
## with this program; if not, write to the Free Software Foundation, Inc.,
## 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import numpy as np
import sys
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.constants import c

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
    if tplname.endswith('.model') or tplname.endswith('.fits'):
        # echelle template
        pixel, wave, spec, err, flag_pixel, bjd, berv = Spectrum(tplname, order=order, targ=targ)
        wave *= 1 + (berv*u.km/u.s/c).to_value('')   # *model already barycentric corrected (?)

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


