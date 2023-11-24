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

import cpl
from cpl.core import Table
from cpl.core import PropertyList, Property


# see https://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py

path = sys.path[0] + '/lib/CRIRES/'

location = crires = EarthLocation.from_geodetic(lat=-24.6268*u.deg,
                                                lon=-70.4045*u.deg,
                                                height=2648*u.m)

oset = '1:19'

ip_guess = {'s': 1.5}

def Spectrum(filename='', order=None, targ=None):
    hdr = PropertyList.load(filename, 0)
    dateobs = hdr["DATE-OBS"].value
    ra = hdr["RA"].value
    de = hdr["DEC"].value

    ndit = hdr["ESO DET NDIT"].value
    nods = hdr["ESO PRO DATANCOM"].value   # Number of combined frames 

    setting = hdr["ESO INS WLEN ID"].value

    hdr = PropertyList.load(filename, 1)
    exptime = hdr["EXPTIME"].value
    naxis = hdr["NAXIS2"].value

    exptime = (exptime*nods*ndit)/2.

    targdrs = SkyCoord(ra=ra*u.deg, dec=de*u.deg)
    if not targ: targ = targdrs
    midtime = Time(dateobs, format='isot', scale='utc') + exptime * u.s
    berv = targ.radial_velocity_correction(obstime=midtime, location=crires)
    berv = berv.to(u.km/u.s).value
    bjd = midtime.tdb

    order_drs, detector = divmod(order-1, 3)
    order_drs = 7 - order_drs		# order number (CRIRES+ definition)
    detector += 1			# detector number (1,2,3)

    tbl = Table.load(filename, detector)
    spec = np.array([tbl["0"+str(order_drs)+"_01_SPEC", i] for i in range(naxis)])
    err = np.array([tbl["0"+str(order_drs)+"_01_ERR", i] for i in range(naxis)])
    pixel = np.arange(spec.size)   

    if str(setting) in ('K2148', 'K2166', 'K2192'):
	# using an own wavelength solution instead of the one created by DRS
        file_wls = np.genfromtxt(path+'wavesolution_own/wave_solution_'+str(setting)+'.dat', dtype=None, names=True).view(np.recarray)
        coeff_wls = [file_wls.b1[order-1], file_wls.b2[order-1], file_wls.b3[order-1]]
        wave = np.poly1d(coeff_wls[::-1])(pixel)

        # using blaze function
        blaze = np.load(path+str(setting)+'_blaze_own.npy')[order]
        spec /= blaze

    else:
        wave = np.array([tbl["0"+str(order_drs)+"_01_WL", i] for i in range(naxis)])
        wave *= 10

    flag_pixel = 1 * np.isnan(spec)		# bad pixel map

    return pixel, wave, spec, err, flag_pixel, bjd, berv


def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''

    if tplname.endswith('_tpl.fits'):
        # tpl created with viper
        hdr = PropertyList.load(tplname, 1)
        naxis = hdr["NAXIS2"].value

        order_drs, detector = divmod(order-1, 3)
        order_drs = 7 - order_drs		# order number (CRIRES+ definition)
        detector += 1			# detector number (1,2,3)

        tbl = Table.load(tplname, detector)
        spec = np.array([tbl["0"+str(order_drs)+"_01_SPEC", i] for i in range(naxis)])
       # err = np.array([tbl["0"+str(order_drs)+"_01_ERR", i] for i in range(naxis)])
        wave = np.array([tbl["0"+str(order_drs)+"_01_WL", i] for i in range(naxis)])
        wave *= 10
    else:
        pixel, wave, spec, err, flag_pixel, bjd, berv = Spectrum(tplname, order=order, targ=targ)
        wave *= 1 + (berv*u.km/u.s/c).to_value('')

    return wave, spec


def FTS(ftsname='lib/CRIRES/FTS/CRp_SGC2_FTStmpl-HR0p007-WN3000-5000_Kband.dat', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)


def write_fits(wtpl_all, tpl_all, e_all, list_files, file_out):

    file_in = list_files[0]

    # copy header from first fits file
    hdr = PropertyList.load(file_in, 0)

    if len(list_files) > 1:
        # delete parts that vary for all observations
        PropertyList.del_regexp(hdr, "DATE-OBS", False)
        PropertyList.del_regexp(hdr, "UTC", False)
        PropertyList.del_regexp(hdr, "LST", False)
        PropertyList.del_regexp(hdr, "ARCFILE", False)
        PropertyList.del_regexp(hdr, "ESO INS SENS*", False)
        PropertyList.del_regexp(hdr, "ESO INS TEMP*", False)
        PropertyList.del_regexp(hdr, "ESO INS1*", False)
        PropertyList.del_regexp(hdr, "ESO DET*", False)
        PropertyList.del_regexp(hdr, "ESO OBS*", False)
        PropertyList.del_regexp(hdr, "ESO TPL*", False)
        PropertyList.del_regexp(hdr, "ESO TEL*", False)
        PropertyList.del_regexp(hdr, "ESO OCS MTRLGY*", False)
        PropertyList.del_regexp(hdr, "ESO ADA*", False)
        PropertyList.del_regexp(hdr, "ESO AOS*", False)
        PropertyList.del_regexp(hdr, "ESO SEQ*", False)
        PropertyList.del_regexp(hdr, "ESO PRO DATANCOM", False)
        PropertyList.del_regexp(hdr, "ESO PRO REC1 PARAM*", False)
        PropertyList.del_regexp(hdr, "ESO PRO REC1 RAW*", False)

    # save raw file informations in FITS header    
    hdr.append(Property('ESO PRO REC2 ID', 'viper_create_tpl', 'Pipeline recipe'))

    for i in range(0, len(list_files), 1):
        pathi, filei = os.path.split(list_files[len(list_files)-i-1])
        hdr.append(Property('ESO PRO REC2 RAW'+str(len(list_files)-i)+' NAME', filei, 'File name'))

    hdr.append(Property('ESO PRO DATANCOM', len(list_files), 'Number of combined frames'))

    for detector in (1, 2, 3):
        # data spread over 3 detectors, each having 6 orders

        # Update headers of the single detectors
        hdro = PropertyList.load(file_in, detector)
        hdro["EXPTIME"].value = 0
        PropertyList.del_regexp(hdro, "ESO *", False)

        tbl = Table.load(file_in, detector)

        for odrs in range(2, 9, 1):    
            o = (7-odrs)*3 + detector
            if o in list(tpl_all.keys()):
                tbl["0"+str(odrs)+"_01_WL"] = wtpl_all[o] / 10.		# wavelength	
                tbl["0"+str(odrs)+"_01_SPEC"] = tpl_all[o]		# data
                tbl["0"+str(odrs)+"_01_ERR"] = e_all[o]			# errors
            else:
               # writing ones for non processed orders
               # tbl["0"+str(odrs)+"_01_WL"] = wtpl_all[o]
                tbl["0"+str(odrs)+"_01_SPEC"] = np.ones(2048)
                tbl["0"+str(odrs)+"_01_ERR"] = np.nan * np.ones(2048)

        if detector == 1:
            Table.save(tbl, hdr, hdro, file_out+'_tpl.fits', cpl.core.io.CREATE)
        else:     
            Table.save(tbl, None, hdro, file_out+'_tpl.fits', cpl.core.io.EXTEND)  

