#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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


oset = '0:56'    # Relative orders that are analyzed (Those are currently ALL AVAILABLE ORDERS)
iset='200:1800'    # Keep the pixels between X and Y ; Do not analyze the ones outside of that range -- Chosen arbitrarily and can be changed by the user
fix = 'wave'    # CARMENES is stabilized : need to -fix wave

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

    ## CARM_NIR has 28 orders but 2 detectors so data is split in half ; data from detectors is merged for obsfiles but not for tpl
    # reshape obsfiles data to match tpl shape : 28 orders * 2 detectors = 56 orders
    dim = spec.shape    # SERVAL tpl shape is (56, 1699) ; obsfile shape is (28, 4080)
    dim = (56, int(dim[1]/(56/dim[0])))    # if tpl : keep same shape ; if obs : change shape to match tpl
    spec = spec.reshape(dim)
    wavelen = wavelen.reshape(dim)
    err_spec = err_spec.reshape(dim)
    
    '''
    good orders
    [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 28, 29, 31, 46, 48, 50, 52]
    '''
    
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
    read_tpl(tplname) reads wavelen and flux from template ; applies necessary corrections depending on type of file (barycentric correc, vacuum correc, ...)
    '''
    wavelen, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ)
    
    # no tpl if tellurics are too strong (= wave, spec set to 0 or nan) ; read_tpl() returns exp(wavelen)
    tpl_no_data = (wavelen==1).all() or np.isnan(wavelen).all()
    if tpl_no_data:
        print(f'\x1b[0;31;40m \tOrder {order} has no template and will not be analyzed. \x1b[0m')
    
    return wavelen, spec



# If using a cell (default = No cell)
def FTS(ftsname='None', dv=100):
    '''FTSFits() reads FTS of the cell and obtains wavenumber and flux (f) from data headers
    Converts wavenumbers [cm] to wavelengths (w) [angstrom], also inverts w and f arrays ( [::-1] ) so that wavelengths are in ascending order
    resample() returns w, f, uj= np.arange(ln(w)) with step dv/c, iod_j=flux interpolated from uj
    '''

    return resample(*FTSfits(ftsname), dv=dv)



# If we want to create a template
def write_fits(wtpl_all, tpl_all, e_all, list_files, file_out):
    
    # Get data from header of the first fits file
    file_in = list_files[0]
    hdu = fits.open(file_in, ignore_blank=True)
    f = hdu[0].data
    
    # Write template data to the file
    for o in range(1, len(f), 1): 
        if o in tpl_all:
            f[o] = tpl_all[o]
        else:
            f[o] = np.ones(len(f[o]))

    hdu.writeto(file_out+'_tpl.model', overwrite=True)


