#! /usr/bin/env python3

## viper - Velocity and IP Estimator
## Copyright (C) Mathias Zechmeister and Jana Koehler
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along
## with this program; if not, write to the Free Software Foundation, Inc.,
## 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import numpy as np

from astropy.io import fits

c = 3e8   # [m/s] speed of light


def FTSfits(ftsname):

    if ftsname.endswith(".dat"):
        data = np.loadtxt(ftsname)
        w = data[:,0]
        f = data[:,1]
        f = f[::-1]
        w = 1e8 / w[::-1]
    else:
        hdu = fits.open(ftsname, ignore_blank=True, output_verify='silentfix')[0]

        f = hdu.data[::-1]
        h = hdu.header
        try:
            w = h['CRVAL1'] + h['CDELT1'] * (np.arange(f.size) + 1. - h['CRPIX1'])
        except:
            # for OES
            w = h['CRVAL1'] + h['CDELT1'] * (np.arange(f.size) + 1.)
        w = 1e8 / w[::-1]   # convert wavenumbers to wavelength [angstrom]

    return w, f


def resample(w, f, dv=100):
    '''
    dv: Sampling step for uniform log(lambda) [m/s]
    '''
    # define a supersampled log(wavelength) space with knot index j
    u = np.log(w)
    uj = np.arange(u[0], u[-1], dv/c)
    iod_j = np.interp(uj, u, f)

    return w, f, uj, iod_j


