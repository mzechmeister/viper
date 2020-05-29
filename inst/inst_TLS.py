from astropy.io import fits
from readmultispec import readmultispec
# see ttps://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py

filename = '../data/TLS/BETA_GEM.fits'
hdu = fits.open(filename)[0]
f = hdu.data[33]
#w_i = 6128.8833940969 + 0.05453566108124*np.arange(f_i.size)  # guess
#w_i, f_i = np.loadtxt('data/BETA_GEM.order34.dat').T

gg = readmultispec(filename, reform=True, quiet=True)
      # "".join(self.header['WAT2_0*'].values()).split("spec")
w = gg['wavelen'] #[orders] # air or vac?

print(w[33])
