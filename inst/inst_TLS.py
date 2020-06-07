import numpy as np
from astropy.io import fits
from .readmultispec import readmultispec
from .airtovac import airtovac

from .FTS_resample import resample, FTSfits

# see https://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py

def Spectrum(filename='data/TLS/other/BETA_GEM.fits', o=None):
    hdu = fits.open(filename)[0]
    f = hdu.data
    gg = readmultispec(filename, reform=True, quiet=True)
    w = gg['wavelen']
    w = airtovac(w)
    if o is not None:
         w, f = w[o], f[o]

    b = 1 * np.isnan(f) # bad pixel map
    b[f>1.5] |= 2 # large flux
    b[(5300<w) & (w<5343)] |= 4  # only for HARPS s1d template (this order misses)

    return w, f, b

def Tpl(tplname, o=None):
    if tplname.endswith('.model'):
        # echelle template
        from inst.inst_TLS import Spectrum
        w, f, b = Spectrum(tplname)
        if o is not None:
            w, f = w[o], f[o]
    if tplname.endswith('_s1d_A.fits'):
        hdu = fits.open(tplname)[0]
        f = hdu.data
        h = hdu.header
        w = h['CRVAL1'] +  h['CDELT1'] * (1. + np.arange(f.size) - h['CRPIX1'])
        w = airtovac(w)
    else:
        # long 1d template
        hdu = fits.open(tplname)
        w = hdu[1].data.field('Arg')
        f = hdu[1].data.field('Fun')

    return w, f


def FTS(ftsname='lib/TLS/FTS/TLS_I2_FTS.fits', dv=100):

    return resample(*FTSfits(ftsname), dv=dv)



