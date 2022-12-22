import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.constants import c

from .readmultispec import readmultispec
from .airtovac import airtovac

from .FTS_resample import resample, FTSfits

# see https://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py

location = oes = EarthLocation.from_geodetic(lat=49.91056*u.deg, lon=14.78361*u.deg, height=528*u.m)

oset = '1:30'

ip_guess = {'s': 300_000/67_000/ (2*np.sqrt(2*np.log(2))) }   # convert FHWM resolution to sigma


def Spectrum(filename='', order=None, targ=None):
    hdu = fits.open(filename, ignore_blank=True)[0]
    hdr = hdu.header

    dateobs = hdr['DATE-OBS']+ 'T' + hdr['UT']
    exptime = hdr['EXPTIME']

    midtime = Time(dateobs, format='isot', scale='utc') + exptime * u.s
    bjd = midtime.tdb

   # targdrs = SkyCoord(ra=ra*u.hour, dec=de*u.deg)
    if not targ: 
        #targ = targdrs
        berv = 0
    else:
        berv = targ.radial_velocity_correction(obstime=midtime, location=oes)
        berv = berv.to(u.km/u.s).value

    spec_obs= hdu.data
    spec_obs/= np.nanmean(spec_obs)
    gg = readmultispec(filename, reform=True, quiet=True)
    wave_obs = gg['wavelen']
    wave_obs = airtovac(wave_obs)
    if order is not None:
         wave_obs, spec_obs= wave_obs[order], spec_obs[order]

    pixel = np.arange(spec_obs.size) 
    e_obs = np.ones(spec_obs.size)*0.1
    flag_pixel = 1 * np.isnan(spec_obs) # bad pixel map
 #   b[f>1.5] |= 4 # large flux

    return pixel, wave_obs, spec_obs, e_obs, flag_pixel, bjd, berv


def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    if tplname.endswith('_s1d_A.fits'):
        hdu = fits.open(tplname)[0]
        spec_tpl= hdu.data
        h = hdu.header
        wave_tpl = h['CRVAL1'] +  h['CDELT1'] * (1. + np.arange(spec_tpl.size) - h['CRPIX1'])
        wave_tpl = airtovac(wave_tpl)
    elif tplname.endswith('1d.fits'):
        hdu = fits.open(tplname, ignore_blank=True)[0]
        hdr = hdu.header
        dateobs = hdr['DATE']
        exptime = hdr['EXPTIME']
        midtime = Time(dateobs, format='isot', scale='utc') + exptime * u.s
        bjd = midtime.tdb

        if not targ: 
            #targ = targdrs
            berv = 0
        else:
            berv = targ.radial_velocity_correction(obstime=midtime, location=oes)
            berv = berv.to(u.km/u.s).value

        spec_tpl = hdu.data
        spec_tpl /= np.nanmean(spec_tpl)
        gg = readmultispec(tplname, reform=True, quiet=True)
        wave_tpl = gg['wavelen']
        wave_tpl = airtovac(wave_tpl)
        wave_tpl *= 1 + (berv*u.km/u.s/c).to_value('')
    else:
        pixel, wave_tpl, spec_tpl, e_tpl, flag_pixel, bjd, berv = Spectrum(tplname, order=order, targ=targ)
        wave_tpl *= 1 + (berv*u.km/u.s/c).to_value('')

    return wave_tpl, spec_tpl


def FTS(ftsname='lib/oes.fits', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)


def Tell(molec):
      modelfile = 'lib/atmos/'+str(molec)+'.fits'
      hdu = fits.open(modelfile, ignore_blank=True)
      atm_model = hdu[1].data
      wave_atm = atm_model.field(0).astype(np.float64)
      spec_atm = atm_model.field(1).astype(np.float64)
    
      return wave_atm, spec_atm


def write_fits(wtpl_all, tpl_all, e_all, list_files, file_out):

    file_in = list_files[0]

    # copy header from first fits file 
    hdu = fits.open(file_in, ignore_blank=True)[0]
    f = hdu.data

    # write the template data to the file
    for order in range(1,49,1): 
        if order in tpl_all:
            f[order] = tpl_all[order]
        else:
            f[order] = np.ones(len(f[order]))

    hdu.writeto(file_out+'.model', overwrite=True)  
  #  hdu.close()  

