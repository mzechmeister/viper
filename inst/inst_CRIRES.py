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
   
    ra = hdr['RA']
    de = hdr['DEC']
    hdr = hdu[1].header
    exptime = hdr['EXPTIME']

    print('---dateobs, exptime  :', dateobs, exptime)

    targdrs = SkyCoord(ra=ra*u.hour, dec=de*u.deg)
    if not targ: targ = targdrs
    midtime = Time(dateobs, format='isot', scale='utc') + exptime * u.s
    berv = targ.radial_velocity_correction(obstime=midtime, location=crires)
    berv = berv.to(u.km/u.s).value
    bjd = midtime.tdb

    oi = 5-int((o-1)/3)  
    d = np.mod(o,3) 
    if d == 0:
        d=3

    w = (hdu[d].data.field(3*oi+2))*10
    e = err = hdu[d].data.field(3*oi+1)
    f = hdu[d].data.field(3*oi)
    x = np.arange(f.size) 

    if (o >= 10) and (o <= 16):
       # using an own wavelength solution as the one given by CRIRES+ pipeline is wrong
        bw = (np.load('wave_solution_cell_best.npy'))[o-1]
        w = np.poly1d(bw[::-1])(x)
    else:
       # still a better wavelength solution for other orders needed
       # should be done with more data in next commisioning run (May 2021)
        obsna = 'data/CRIRES/210219_Ross619/withSGC/cr2res_obs_nodding_extracted_combined.fits'
        hdu = fits.open(filename, ignore_blank=True)
        w = (hdu[d].data.field(3*oi+2))*10
        shift_cell = [111.6,110.2,80.76,88.12,70,120,138,82.1,81.02,80.98,51.67, 81.78,102.48, 66.35,98.11,92.08,64.73,110.06]       
        w /= (1+shift_cell[o-1]/3e5)
 

    # correction of instrumental shift
    # applying a offset correction for wavelength solution
    # is it needed? is it right?
 #   sss = [ 2.59857143-0.17,  2.80714286-0.22,  3.21571429-0.2,  3.03285714+0.58,  3.64428571-0.26,  3.45857143, 1.90571429+0.37]
 #   sss = [1.944675, 2.206255, 2.5936, 3.1737325, 3.14364, 3.1717725, 2.3365325]
  #  w /= (1+(sss[o-10]+0.2)/3e5)

    print("---berv, wlmin, wlmax:", berv, w[0],w[-1])
    
    fwhm = (hdu[d]).header['ESO QC SLITFWHM'+str(oi+2)]
    snr = (hdu[d]).header['ESO QC SNR'+str(oi+2)]
    slit = (hdu[0]).header['ESO INS SLIT1 WID']

    print("---fwhm, snr, slit   :", fwhm, snr, slit)

#    w = airtovac(w)	# still needed with own ws?
    
    b = 1 * np.isnan(f) # bad pixel map
 #   b[f>15000] |= 2 # large flux

    return x, w, f, b, bjd, berv

def Tpl(tplname, o=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    if tplname.endswith('.model'):
        # echelle template
        x, w, f, b, bjd, berv = Spectrum(tplname, o=o, targ=targ)
        w *= 1 + (berv*u.km/u.s/c).to_value('')   # *model already barycentric corrected (?)
    else:
        # long 1d template
        x, w, f, b, bjd, berv = Spectrum(tplname, o=o, targ=targ)   
        w *= 1 + (berv*u.km/u.s/c).to_value('')

    return w, f


def FTS(ftsname='lib/CRIRES/FTS/CRp_SGC2_FTStmpl-HR0p007-WN3000-5000_Kband.dat', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)
