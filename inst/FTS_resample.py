import numpy as np
from scipy import interpolate

from astropy.io import fits

c = 3e8   # [m/s] speed of light


def FTSfits(ftsname):
    hdu = fits.open(ftsname)[0]

    f = hdu.data[::-1]
    h = hdu.header
    w = h['CRVAL1'] + h['CDELT1'] * (np.arange(f.size) + 1. - h['CRPIX1'])
    w = 1e8 / w[::-1]   # convert wavenumbers to wavelength [angstrom]

    return w, f


def resample(w, f, dv=100):
    '''
    dv: Sampling step for uniform log(lambda) [m/s]
    '''
    # define a supersampled log(wavelength) space with knot index j
    xj = np.arange(np.log(w[0]), np.log(w[-1]), dv/c)
    iod_j = interpolate.interp1d(np.log(w), f)(xj)

    return w, f, xj, iod_j


