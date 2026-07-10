#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:10:29 2026

This file was written by Naïm Teriitehau-Martin. For any issues, please contact : naimteriitehau@gmail.com 
"""

import numpy as np
import os
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.constants import c

from inst.template import read_tpl
from inst.FTS_resample import resample, FTSfits 

'''
If creating or using a template with viper : oset = '0:56'
If using a template from SERVAL : oset = '3:56'
'''
oset = '0:56'    # Relative orders that are analyzed -- Trying to ignore bad orders
iset='200:1800'    # Keep the pixels between X and Y ; Do not analyze the ones outside of that range -- Chosen arbitrarily and can be changed by the user
fix = 'wave ip'    # CARMENES is stabilized : need to -fix wave ; fixing the IP gives better results

# Convert FWHM resolution to sigma
ip_guess = {'s' : c.to_value(u.km/u.s)/(80_400*2*np.sqrt(2*np.log(2)))}  # speed of light [km/s] divided by (spectrograph_resolution x factor) -- Assuming Gaussian IP for factor ; FWHM is defined in terms of speed dv [km/s] (which is ~= delta(ln(wavelen))) 

# Location of CARMENES Spectrograph -- Obtained from the data headers (hdr keys 'HIERARCH CAHA TEL GEOLAT' and 'HIERARCH CAHA TEL GEOLON' and 'HIERARCH CAHA TEL GEOELEV')
location = carmenes = EarthLocation.from_geodetic(lat=37.2236*u.deg, lon=-2.54625*u.deg, height=2168.*u.m)


def Spectrum(filename='', order=None, targ=None):
    hdu = fits.open(filename, ignore_blank=True)
    hdr = hdu[0].header
    
    dateobs = hdr.get('DATE-OBS')   # UTC datetime at observation start (ISOT format)
    exptime_tmean = hdr.get('HIERARCH CARACAL TMEAN') * u.s  # Flux-weighted midpoint of exposure

    # If target not specified while calling Spectrum() : Define target from data file header info
    ra = hdr.get('RA', np.nan)     # RightAscension of target [degrees]
    dec = hdr.get('DEC', np.nan)    # Declination of target [degrees]
    dec = (dec+90) % 180 - 90   #Declination needs to be between -90 and 90 deg, but some SERVAL tpl headers give values way bigger than that so we fix it
    targdrs = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    if not targ: targ = targdrs
    
    
    # Apply barycentric RV correction at exposure flux-weighted midpoint
    midtime = Time(dateobs, format='isot', scale='utc') + exptime_tmean
        
    berv = targ.radial_velocity_correction(obstime=midtime, location=carmenes)  #Barycentric Earth RV : This is the RV correction to apply
    berv = berv.to_value(u.km/u.s)
    bjd = midtime.tdb    # Convert midtime scale from utc to Barycentric Dynamical Time
    
    
    # Read file data to obtain spec, wavelen, err
    spec = hdu['SPEC'].data    # Flux for spectrum
    wavelen = hdu['WAVE'].data    # Vacuum wavelength for each pixel [angstrom]
    err_spec = hdu['SIG'].data    # Error estimate for flux
    
    
    ## CARM_NIR has 28 orders but 2 detectors so data is split in half. We need to split each order in half.
    # reshape arrays : 28 orders * 2 detectors = 56 orders
    dim = spec.shape
    dim = (56, int(dim[1]/(56/dim[0])))    # Split into 56 orders of equal size (e.g. order 0 is split into order 0 and order 1)
    spec = spec.reshape(dim)
    wavelen = wavelen.reshape(dim)
    err_spec = err_spec.reshape(dim)
    
        
    # If specific order is selected, only keep that order
    if order is not None:
        wavelen, spec, err_spec = wavelen[order], spec[order], err_spec[order]
    
    # Build pixel array and bad pixel map
    pixel = np.arange(spec.size)    # As many pixels as there are values in spec (Including bad pixels)
    flag_pixel = 1 * np.isnan(spec)  # Bad pixels (where spec = nan) have a value of 1
    
    return pixel, wavelen, spec, err_spec, flag_pixel, bjd, berv



def Tpl(tplname, order=None, targ=None):
    '''
    Tpl should return barycentric corrected wavelengths.
    read_tpl(tplname) reads wavelen and flux from template file ; applies necessary corrections depending on type of file (barycentric correc, vacuum correc, ...)
    '''
    wavelen, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ)
    
    # if no tpl data (e.g. if tellurics are too strong) : wave, spec are set to 0 or nan ; read_tpl() returns exp(wavelen)
    tpl_no_data = (wavelen==1).all() or np.isnan(wavelen).all()
    if tpl_no_data:
        #wavelen = np.full_like(wavelen, np.nan)
        #spec = np.full_like(wavelen, np.nan)
        print(f'\x1b[0;31;40m \tOrder {order} has no template and will not be analyzed. \x1b[0m')
    
    return wavelen, spec



# If using a cell (default = No cell) -- CARMENES does not use a cell so we should not need this
def FTS(ftsname='None', dv=100):
    '''FTSFits() reads FTS of the cell and obtains wavenumber and flux (f) from data headers
    Converts wavenumbers [cm] to wavelengths (w) [angstrom], also inverts w and f arrays ( [::-1] ) so that wavelengths are in ascending order
    resample() returns w, f, uj= np.arange(ln(w)) with step dv/c, iod_j=flux interpolated from uj
    '''

    return resample(*FTSfits(ftsname), dv=dv)




# For template creation
def write_fits(wave_tpl_all, spec_tpl_all, err_tpl_all, list_files, file_out):
    '''
    When creating a template : Orders 17-20 and 36-45 usually give poor results due to very strong tellurics
    '''
    
    # Get data from header of the first fits file
    file_in = list_files[0]
    hdu = fits.open(file_in, ignore_blank=True)
    hdr = hdu[0].header

    # grab shape of data arrays
    spec = hdu['SPEC'].data
    # reshape 28 orders (1 order split between 2 detectors) --> 56 orders
    dim = spec.shape
    dim = (56, int(dim[1]/(56/dim[0])))    # e.g. order N is split into order 2N (from detector1) and order 2N+1 (from detector2)
    
    # set everything to nan by default (for orders that are not in the tpl)
    wavelen = np.full(dim, np.nan)   # Vacuum wavelength for each pixel [angstrom]
    spec = np.full(dim, np.nan)      # Flux for spectrum
    err_spec = np.full(dim, np.nan)  # Error estimate for flux
    
    # write template data to the file
    for o in list(spec_tpl_all.keys()):
        wavelen[o] = wave_tpl_all[o]
        spec[o] = spec_tpl_all[o]
        err_spec[o] = err_tpl_all[o]
    
    hdu['SPEC'].data = spec   # Flux for spectrum
    hdu['WAVE'].data = wavelen    # Vacuum wavelength for each pixel [angstrom]
    hdu['SIG'].data = err_spec    # Error estimate for flux
    


    ### Add header cards
    hdr['DATE-TPL'] = (Time.now().isot , 'UT date of template creation')
    
    for n, obsname in enumerate(list_files):
        hdr[f'HIERARCH VIPER COADD FILE {n}'] = os.path.basename(obsname)    # name of coadded files
        
    hdr['HIERARCH VIPER COADD COMIN'] = (min(spec_tpl_all.keys()), 'minimum coadded order')
    hdr['HIERARCH VIPER COADD COMAX'] = (max(spec_tpl_all.keys()), 'maximum coadded order')
    
    hdr['HIERARCH VIPER COADD NUM'] = (len(list_files), 'number of coadded spectra')
    
    
   # hdu.writeto needs to be fixed because some cards are in the wrong order
    ## --> The fix puts EXTNAME(='SIG','WAVE','SPEC') after PCOUNT and GCOUNT instead of before
    hdu.writeto(file_out+'_tpl.fits', output_verify='silentfix', overwrite=True)
    hdu.close()

