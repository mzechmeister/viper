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

pg = {'s': 300_000/67_000/ (2*np.sqrt(2*np.log(2))) }   # convert FHWM resolution to sigma

def Spectrum(filename='', o=None, targ=None):
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

    f = hdu.data
    f /= np.nanmean(f)
    gg = readmultispec(filename, reform=True, quiet=True)
    w = gg['wavelen']
    w = airtovac(w)
    if o is not None:
         w, f = w[o], f[o]

    x = np.arange(f.size) 
    b = 1 * np.isnan(f) # bad pixel map
 #   b[f>1.5] |= 4 # large flux

    return x, w, f, b, bjd, berv

def Tpl(tplname, o=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    if tplname.endswith('_s1d_A.fits'):
        hdu = fits.open(tplname)[0]
        f = hdu.data
        h = hdu.header
        w = h['CRVAL1'] +  h['CDELT1'] * (1. + np.arange(f.size) - h['CRPIX1'])
        w = airtovac(w)
    else:
        x, w, f, b, bjd, berv = Spectrum(tplname, o=o, targ=targ)
        w *= 1 + (berv*u.km/u.s/c).to_value('')

    return w, f

def FTS(ftsname='lib/oes.fits', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)


