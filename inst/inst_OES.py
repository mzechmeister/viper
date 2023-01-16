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

    spec = hdu.data
    spec /= np.nanmean(spec)
    gg = readmultispec(filename, reform=True, quiet=True)
    wave = gg['wavelen']
    wave = airtovac(wave)
    if order is not None:
         wave, spec= wave[order], spec[order]

    pixel = np.arange(spec.size) 
    err = np.ones(spec.size)*0.1
    flag_pixel = 1 * np.isnan(spec) # bad pixel map
 #   b[f>1.5] |= 4 # large flux

    return pixel, wave, spec, err, flag_pixel, bjd, berv


def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    if tplname.endswith('_s1d_A.fits'):
        hdu = fits.open(tplname)[0]
        spec= hdu.data
        h = hdu.header
        wave = h['CRVAL1'] +  h['CDELT1'] * (1. + np.arange(spec.size) - h['CRPIX1'])
        wave = airtovac(wave)
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

        spec = hdu.data
        spec /= np.nanmean(spec)
        gg = readmultispec(tplname, reform=True, quiet=True)
        wave = gg['wavelen']
        wave = airtovac(wave)
        wave *= 1 + (berv*u.km/u.s/c).to_value('')
    else:
        pixel, wave, spec, err, flag_pixel, bjd, berv = Spectrum(tplname, order=order, targ=targ)
        wave *= 1 + (berv*u.km/u.s/c).to_value('')

    return wave, spec


def FTS(ftsname='lib/oes.fits', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)


def Tell(molec='all'):

      if molec[0] == 'all':
           molec = ['H2O', 'O2']

      wave_atm_all, specs_molec_all = {}, {}
      for mol in range(len(molec)):
          modelfile = 'lib/atmos/'+str(molec[mol])+'.fits'
          hdu = fits.open(modelfile, ignore_blank=True)
          atm_model = hdu[1].data
          w_atm = atm_model.field(0).astype(np.float64)
          f_atm = atm_model.field(1).astype(np.float64)
          wave_atm_all[mol], specs_molec_all[mol] = w_atm, f_atm
    
      return wave_atm_all, specs_molec_all, molec


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

