import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.constants import c

from .readmultispec import readmultispec
from .airtovac import airtovac

from .FTS_resample import resample, FTSfits

import matplotlib.pyplot as plt


# see https://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py

location = crires = EarthLocation.from_geodetic(lat=-24.6268*u.deg, lon=-70.4045*u.deg, height=2648*u.m)

oset = '0:18'

#1.07
# signal: 17394
#slitfwhm = 10.2
pg = {'s': 300_000/17_394/ (2*np.sqrt(2*np.log(2))) }   # convert FHWM resolution to sigma
pg = {'s': 2 }

def Spectrum(filename='', o=None, targ=None):
    hdu = fits.open(filename, ignore_blank=True)
    hdr = hdu[0].header
    dateobs = hdr['DATE-OBS']
   
    ra = hdr.get('RA',np.nan)
    de = hdr.get('DEC',np.nan)
    hdr = hdu[1].header
    exptime = hdr.get('EXPTIME',0) 

    targdrs = SkyCoord(ra=ra*u.deg, dec=de*u.deg)
    if not targ: targ = targdrs
    midtime = Time(dateobs, format='isot', scale='utc') + exptime * u.s
    berv = targ.radial_velocity_correction(obstime=midtime, location=crires)
    berv = berv.to(u.km/u.s).value
    bjd = midtime.tdb

    oi = 5-int((o-1)/3)  
    d = np.mod(o,3) 
    if d == 0:
        d = 3

    e = err = hdu[d].data.field(3*oi+1)
    f = hdu[d].data.field(3*oi)
    x = np.arange(f.size) 

    setting = (hdu[0]).header['ESO INS WLEN ID']
    # using an own wavelength solution as the one given by CRIRES+ pipeline is imprecisely
    if str(setting) in ('K2148','K2166','K2192'):
        B = np.genfromtxt('lib/CRIRES/wavesolution/wave_solution_'+str(setting)+'.dat', dtype=None, names=True).view(np.recarray)
        bw = [B.b1[o-1], B.b2[o-1], B.b3[o-1]]
        w = np.poly1d(bw[::-1])(x)
        # normalization
        bl = np.load('lib/CRIRES/'+str(setting)+'_blaze.npy')[o]
        f /= bl
    else:
        w = (hdu[d].data.field(3*oi+2))*10

    b = 1 * np.isnan(f) # bad pixel map

    # using flag file for noisy regions
    A = np.genfromtxt('lib/CRIRES/flag_file.dat', dtype=None, names=True).view(np.recarray)
    flag_start = A.ok_s
    flag_end = A.ok_e

    b[:flag_start[o-1]] |= 64
    b[flag_end[o-1]:] |= 64

    return x, w, f, b, bjd, berv

def Tpl(tplname, o=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    if tplname.endswith('.npy'):
        w = np.load(tplname)[o,0]
        f = np.load(tplname)[o,1]
    else:
        # long 1d template
        x, w, f, b, bjd, berv = Spectrum(tplname, o=o, targ=targ)   
        w *= 1 + (berv*u.km/u.s/c).to_value('')

    return w, f


def FTS(ftsname='lib/CRIRES/FTS/CRp_SGC2_FTStmpl-HR0p007-WN3000-5000_Kband.dat', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)

def Tell(molec):
      modelfile = 'lib/CRIRES/atmos/stdAtmos_crires_'+str(molec)+'.fits'
      hdu = fits.open(modelfile, ignore_blank=True)
      atm_model = hdu[1].data
      w_atm = atm_model.field(0).astype(np.float64)
      f_atm = atm_model.field(1).astype(np.float64)
    
      return w_atm, f_atm









