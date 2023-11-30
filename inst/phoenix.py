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

import os
from astropy.io import fits

def read(fitsphoe, wmin=3500, wmax=8000):
    '''
    Example:
    --------
    >>> read('lte05100-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')

    The wave file must be in the same directory.

    Phoenix spectra might be used to get absolute RVs, Teff, logg and [Fe/H].
    Also a vsini broadening would be needed (see https://github.com/mzechmeister/serval/blob/a348b4ca77e57b0e9e626f8c9fb147f080cc2418/src/serval.py#L228).
    Of course, the precision with depend on the model (mis-)match of the Phoenix spectra.
    '''
    print('read phoenix', fitsphoe)
    with fits.open(fitsphoe) as hdulist:
        flux = hdulist[0].data
    with fits.open(os.path.join(os.path.dirname(fitsphoe), 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')) as hdulist:
        wave = hdulist[0].data
    imap = (wmin<wave) & (wave<wmax)
    return wave[imap], flux[imap]

