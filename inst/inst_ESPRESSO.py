#! /usr/bin/env python3
# Licensed under a GPLv3 style license - see LICENSE

import numpy as np
import os.path
import sys
from datetime import datetime
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.constants import c

from .readmultispec import readmultispec
from .airtovac import airtovac

from .FTS_resample import resample, FTSfits


# see https://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py


location = espresso = EarthLocation.from_geodetic(
    lat=-24.6276 * u.deg, lon=-70.4051 * u.deg, height=2648 * u.m
)

oset = '40:170'


ip_guess = {'s': 1.5}

def Spectrum(filename='', order=None, targ=None):

    hdu = fits.open(filename, ignore_blank=True)
    hdr = hdu[0].header
    ra = hdr.get('RA', np.nan)
    de = hdr.get('DEC', np.nan)
        
    dateobs = hdr['DATE-OBS']
    berv = hdr['ESO QC BERV']

    spec = hdu['SCIDATA'].data[order]
    wave = hdu['WAVEDATA_VAC_BARY'].data[order]
    err = hdu['ERRDATA'].data[order]
        
    # wavelengths are already berv corrected
    # leads to problems in the telluric corrections for large berv
    wave *= 1-berv/3e5

    pixel = np.arange(spec.size)

  #  targdrs = SkyCoord(ra=ra*u.deg, dec=de*u.deg)
   # if not targ: targ = targdrs
    midtime = Time(dateobs, format='isot', scale='utc') 
 #   berv = targ.radial_velocity_correction(obstime=midtime, location=espresso)
   # berv = 0#berv.to(u.km/u.s).value
    bjd = midtime.tdb

    flag_pixel = 1 * np.isnan(spec)		# bad pixel map
    flag_pixel[spec==1] |= 1
    flag_pixel[spec<=0] |= 1

    return pixel, wave, spec, err, flag_pixel, bjd, berv


def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''

    pixel, wave, spec, err, flag_pixel, bjd, berv = Spectrum(tplname, order=order, targ=targ)
    wave *= 1 + (berv*u.km/u.s/c).to_value('')

    return wave, spec


def FTS(ftsname='', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)


def write_fits(wtpl_all, tpl_all, e_all, list_files, file_out):

    file_in = list_files[0]

    # copy header from first fits file
    hdu = fits.open(file_in, ignore_blank=True)
    hdr = hdu[0].header

    if len(list_files) > 1:
        # delete parts that vary for all observations
        del hdr['DATE-OBS']
        del hdr['UTC']
        del hdr['LST']
        del hdr['ARCFILE']
        del hdr['ESO INS SENS*']
        del hdr['ESO INS TEMP*']
        del hdr['ESO INS1*']
        del hdr['ESO DET*']
        del hdr['ESO OBS*']
        del hdr['ESO TPL*']
        del hdr['ESO TEL*']
        del hdr['ESO OCS MTRLGY*']
        del hdr['ESO ADA*']
        del hdr['ESO AOS*']
        del hdr['ESO SEQ*']
        del hdr['ESO PRO DATANCOM']
        del hdr['ESO PRO REC1 PARAM*']
        del hdr['ESO PRO REC1 RAW*']

        for hdri in hdu:
            hdri.header['EXPTIME'] = 0

    # file creation date
    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%dT%H:%M:%S")
    hdr['DATE'] = dt_string

    # save raw file informations in FITS header
    hdr.set('HIERARCH ESO PRO REC2 ID', 'viper_create_tpl', 'Pipeline recipe', after='ESO PRO REC1 PIPE ID')

    for i in range(0, len(list_files), 1):
        pathi, filei = os.path.split(list_files[len(list_files)-i-1])
        hdr.set('HIERARCH ESO PRO REC2 RAW'+str(len(list_files)-i)+' NAME', filei, 'File name', after='ESO PRO REC2 ID')

    hdr.set('HIERARCH ESO PRO DATANCOM2', len(list_files), 'Number of combined frames', after='ESO PRO REC2 RAW'+str(len(list_files))+' NAME')
    
    spec = hdu['SCIDATA'].data
    wave = hdu['WAVEDATA_VAC_BARY'].data
    err = hdu['ERRDATA'].data

    # write the template data to the file            
    for o in range(0, 170, 1):  

        if o in list(tpl_all.keys()):     
            wave[o] = wtpl_all[o]		# wavelength	
            spec[o] = tpl_all[o]		# data
            err[o] = e_all[o]			# errors   
        else:
           # writing ones for non processed orders
            pix = len(spec[o])
        #    wave[o] = wtpl_all[o]		# wavelength	
            spec[o] = np.ones(pix)		# data
            err[o] = np.nan * np.ones(pix)
            
    hdu.writeto(file_out+'_tpl.fits', overwrite=True)
    hdu.close()
