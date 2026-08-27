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


location = xshooter = EarthLocation.from_geodetic(
    lat=-24.6268 * u.deg, lon=-70.4045 * u.deg, height=2648 * u.m
)

oset = '1:12'

ip_guess = {'s': 300_000/20_000/ (2*np.sqrt(2*np.log(2))) } 

def Spectrum(filename='', order=None, targ=None):

    hdu = fits.open(filename, ignore_blank=True)
    hdr = hdu[0].header
    ra = hdr.get('RA', np.nan)
    de = hdr.get('DEC', np.nan)

    dateobs = hdr['DATE-OBS']
    berv = hdr['ESO QC VRAD BARYCOR']

	# ESO Phase 3 format
	# entire 1D spectra - separation into smaller chunks for viper
    spec = hdu[1].data['FLUX'][0]
    wave = hdu[1].data['WAVE'][0]
    err = hdu[1].data['ERR'][0]
    
    if 1:
        px = 2048	# pixel number of original data ?
        spec = spec[order*px:(order+1)*px]
        wave = wave[order*px:(order+1)*px]*10
        err = err[order*px:(order+1)*px]

 #   wave *= 1-berv/3e5		# wavelengths berv corrected?
    
    wave = airtovac(wave)

    midtime = Time(dateobs, format='isot', scale='utc') 
    bjd = midtime.tdb

    pixel = np.arange(spec.size)
    flag_pixel = 1 * np.isnan(spec)		# bad pixel map
    flag_pixel[spec==1] |= 1
    flag_pixel[spec<=0] |= 1

    return pixel, wave, spec, err, flag_pixel, bjd, berv


def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''

    wave, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ)

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
