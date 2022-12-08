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

location = tls = EarthLocation.from_geodetic(lat=50.980111*u.deg, lon=11.711167*u.deg, height=342*u.m)

oset = '18:30'
iset = '380:1700'
#iset = '0:2048'

pg = {'s': 300_000/67_000/ (2*np.sqrt(2*np.log(2))) }   # convert FHWM resolution to sigma

def Spectrum(filename='data/TLS/other/BETA_GEM.fits', o=None, targ=None):
    hdu = fits.open(filename, ignore_blank=True)[0]
    hdr = hdu.header

    dateobs = hdr['DATE-OBS']
    exptime = hdr.get('EXP_TIME',hdr.get('EXPOSURE'))   # 20211018_guenther_TCEcell_0063.fits EXPOSURE (no exptime)
    ra = hdr.get('RA', np.nan)                          # 20211018_guenther_TCEcell_0184.fits no RA
    de = hdr.get('DEC', np.nan)

    targdrs = SkyCoord(ra=ra*u.hour, dec=de*u.deg)
    if not targ: targ = targdrs
    midtime = Time(dateobs, format='isot', scale='utc') + exptime * u.s
    berv = targ.radial_velocity_correction(obstime=midtime, location=tls)
    berv = berv.to(u.km/u.s).value
    bjd = midtime.tdb

    f = hdu.data
    gg = readmultispec(filename, reform=True, quiet=True)
    w = gg['wavelen']
    w = airtovac(w)
    if o is not None:
         w, f = w[o], f[o]

    x = np.arange(f.size) 
    e = np.zeros(f.size)+0.1
    b = 1 * np.isnan(f) # bad pixel map
    #b[f>1.] |= 4   # large flux, only for normalised spectra, use kapsig instead
    b[(5300<w) & (w<5343)] |= 256  # only for HARPS s1d template (this order misses)
    # TLS spectra have a kink in continuum  at about 1700
    # Also the deconv could have a bad wavelength solution.

    return x, w, f,e, b, bjd, berv

def Tpl(tplname, o=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    if tplname.endswith('.model'):
        # echelle template
        x, w, f, e, b, bjd, berv = Spectrum(tplname, o=o, targ=targ)
        w *= 1 + (berv*u.km/u.s/c).to_value('')   # *model already barycentric corrected (?)
    elif tplname.endswith('_s1d_A.fits') or  tplname.endswith('.tpl.s1d.fits'):
        hdu = fits.open(tplname)[0]
        f = hdu.data
        h = hdu.header
        w = h['CRVAL1'] +  h['CDELT1'] * (1. + np.arange(f.size) - h['CRPIX1'])
        if tplname.endswith('_s1d_A.fits'):
            w = airtovac(w)
        else:
            w = np.exp(w)
    elif tplname.endswith('PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'):
        from . import phoenix
        w, f = phoenix.read(tplname)
    else:
        # long 1d template
        hdu = fits.open(tplname)
        w = hdu[1].data.field('Arg')
        f = hdu[1].data.field('Fun')

    return w, f


def FTS(ftsname='lib/TLS/FTS/TLS_I2_FTS.fits', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)

def Tell(molec):
      modelfile = 'lib/atmos/'+str(molec)+'.fits'
      hdu = fits.open(modelfile, ignore_blank=True)
      atm_model = hdu[1].data
      w_atm = atm_model.field(0).astype(np.float64)
      f_atm = atm_model.field(1).astype(np.float64)
    
      return w_atm, f_atm

def write_fits(wtpl_all, tpl_all, e_all, list_files, file_out):

    file_in = list_files[0]

    # copy header from first fits file 
    hdu = fits.open(file_in, ignore_blank=True)[0]
    f = hdu.data

    # write the template data to the file
    for o in range(1,49,1): 
        if o in tpl_all:
            f[o] = tpl_all[o]
        else:
            f[o] = np.ones(2048)

    hdu.writeto(file_out+'.model', overwrite=True) 


