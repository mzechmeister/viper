import numpy as np
from astropy.io import fits
from scipy import interpolate


c = 3e5   # [km/s] speed of light
def FTS(ftsname='lib/TLS/FTS/TLS_I2_FTS.fits'):
    hdu = fits.open(ftsname)[0]

    f = hdu.data[::-1]
    h = hdu.header
    w = h['CRVAL1'] + h['CDELT1'] * (np.arange(f.size) + 1. - h['CRPIX1'])
    w = 1e8 / w[::-1]   # convert wavenumbers to wavelength [angstrom]

    xj = np.arange(np.log(w[0])+100/c, np.log(w[-1])-100/c, 0.1/c)# reduce range by 100 km/s
    iod_j = interpolate.interp1d(np.log(w), f)(xj)

    
    dx = xj[1] - xj[0]  
    print("sampling [km/s]:", dx*c)
    return w,f,xj, iod_j



