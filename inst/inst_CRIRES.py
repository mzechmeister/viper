import numpy as np
import os.path
from datetime import datetime
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

oset = '1:19'

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

    oi, d = divmod(o-1, 3)
    oi = 5 - oi	# order number (CRIRES+ definition)
    d += 1		# detector number (1,2,3)

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

    return x, w, f, e, b, bjd, berv

def Tpl(tplname, o=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''

    if tplname.endswith('_tpl.fits'):
        # tpl created with viper
        hdu = fits.open(tplname, ignore_blank=True)
        hdr = hdu[0].header

        oi, d = divmod(o-1, 3)
        oi = 5 - oi	# order number (CRIRES+ definition)
        d += 1		# detector number (1,2,3)

        e = err = hdu[d].data.field(3*oi+1)
        f = hdu[d].data.field(3*oi)
        x = np.arange(f.size) 
        w = (hdu[d].data.field(3*oi+2))
    elif tplname.endswith('.npy'):
        w = np.load(tplname)[o,0]
        f = np.load(tplname)[o,1]
    else:
        # long 1d template
        x, w, f, e, b, bjd, berv = Spectrum(tplname, o=o, targ=targ)   
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
      # add wavelength shift, as synthetic telluric spectra are laboratory wavelengths
      # shift was determined empirical on several observations
      w_atm *= (1 + (-0.249/3e5))
    
      return w_atm, f_atm


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
    hdr.set('ESO PRO REC2 ID', 'viper_create_tpl', 'Pipeline recipe', after='ESO PRO REC1 PIPE ID')

    for i in range(0, len(list_files), 1):
        pathi, filei = os.path.split(list_files[len(list_files)-i-1])
        hdr.set('ESO PRO REC2 RAW'+str(len(list_files)-i)+' NAME', filei, 'File name', after='ESO PRO REC2 ID')

    hdr.set('ESO PRO DATANCOM', len(list_files), 'Number of combined frames', after='ESO PRO REC2 RAW'+str(len(list_files))+' NAME')

    # write the template data to the file
    for o in range(1,19,1): 
        # data spread over 3 detectors, each having 6 orders
        oi, d = divmod(o-1, 3)
        oi = 5 - oi	# order number (CRIRES+ definition)
        d += 1		# detector number (1,2,3)

        data = hdu[d].data
        cols = hdu[d].columns

        if o in list(tpl_all.keys()):
            data[str(cols.names[3*oi])] = tpl_all[o]			# data
            data[str(cols.names[3*oi+1])] = e_all[o]			# errors
            data[str(cols.names[3*oi+2])] = wtpl_all[o]			# wavelength
        else:
            # writing zeros for non processed orders
            data[str(cols.names[3*oi])] = np.ones(2048)
            data[str(cols.names[3*oi+1])] = np.nan * np.ones(2048)
            data[str(cols.names[3*oi+2])] = (data.field(3*oi+2))*10	 # [Angstrom]

    hdu.writeto(file_out+'_tpl.fits', overwrite=True)  
    hdu.close()  







