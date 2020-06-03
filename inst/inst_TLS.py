from astropy.io import fits
from .readmultispec import readmultispec
from .airtovac import airtovac
# see ttps://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py

def Spectrum(filename = 'data/TLS/other/BETA_GEM.fits'):
   hdu = fits.open(filename)[0]
   f = hdu.data
   gg = readmultispec(filename, reform=True, quiet=True)
   w = gg['wavelen'] #[orders] # [air]
   w = airtovac(w)
   return w, f

