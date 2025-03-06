#! /usr/bin/env python3
# Licensed under a GPLv3 style license - see LICENSE

import numpy as np
import sys
import os
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

from .template import read_tpl
from .airtovac import airtovac

from .FTS_resample import resample, FTSfits
from pause import pause

location = lasilla = EarthLocation.of_site('lasilla')

oset = '0:1'

ip_guess = {'s': 300_000/200_000/ (2*np.sqrt(2*np.log(2))) }   # convert FHWM resolution to sigma

def Spectrum(filename, order=None, targ=None, chksize=4000):
    with open(filename) as myfile:
        hdr = [next(myfile) for x in range(21)]
    # PX#   WAVELENGTH          FLUX           ERROR         MASK (0/1/6)
    x, w, f, e_f, m = np.genfromtxt(filename, skip_header=21).T
    w = airtovac(w)
    if 1:
        # stitching
        ax = 0*w + 1    # default pixel size
        ax[576] = 0.963
        ax[4+512*1] = 0.986
        ax[4+512*2-1] = 0.981
        ax[4+512*2] = 0.961
        ax[4+512*3-1] = 0.993
        ax[4+512*4] = 0.987
        ax[4+512*5] = 0.992
        ax[4+512*6-1] = 0.986
        ax[4+512*7] = 0.971
        x = np.cumsum(ax) - 1

    if order is not None:
        order = slice(order*chksize, (order+1)*chksize)
        x, w, f, e_f, m = x[order], w[order], f[order], e_f[order], m[order]

    b = 1 * np.isnan(f) # bad pixel map
    #b[f>1.5] |= 2 # large flux
    b[2120:2310] |= 4  # grating ghost on spectrum, CES.2000-08-13T073047.811, stationary?

    dateobs = hdr[2].split()[-1]
    exptime = float(hdr[4].split()[-1])
    if ':' not in dateobs:  # e.g. HR2667_1999-12-30T062739.773.dat
         dateobs = hdr[0].split()[-1]
         exptime = float(hdr[2].split()[-1])

    ra = '03:17:46.1632605674'
    de = '-62:34:31.154247481'
    ra = '03:18:12.8185412558'
    de = '-62:30:22.917300282'
    pmra = 1331.151
    pmde = 648.523
    #from pause import pause; pause()
    #SkyCoord.from_name('M31', frame='icrs')
    # sc = SkyCoord(ra=ra, dec=de, unit=(u.hourangle, u.deg), pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmde*u.mas/u.yr)
    midtime = Time(dateobs, format='isot', scale='utc') + exptime/2 * u.s
    berv = targ.radial_velocity_correction(obstime=midtime, location=lasilla)
    berv = berv.to(u.km/u.s).value
    bjd = midtime.tdb
    return x, w, f, e_f, b, bjd, berv

def Tpl(tplname, order=None, targ=None):
    if tplname.endswith('.dat'):
        # echelle template
        x, w, f, e, b, bjd, berv = Spectrum(tplname, targ=targ)
        w *= 1 + berv/3e5
    else:
        wave, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ)

    return w, f


#def FTS(ftsname='lib/CES/iodine_50_wn.fits', dv=100):
def FTS(ftsname='lib/CES/master_iodlund.dat', dv=100):
    if ftsname.endswith('iodlund.dat'):
        w = airtovac(1e8 / np.fromfile(ftsname)[-2::-2])
        f = np.fromfile(ftsname)[::-2]
        return resample(w, f, dv=dv)
    else:
        print('FTS', ftsname)
        return resample(*FTSfits(ftsname), dv=dv)


