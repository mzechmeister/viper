from astropy.io import fits
import readmultispec
# see ttps://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py

filename = '../TLS/data/BETA_GEM.fits'
hdu = fits.open(filename)[0]
f_i = hdu.data[33]
#w_i = 6128.8833940969 + 0.05453566108124*np.arange(f_i.size)  # guess
#w_i, f_i = np.loadtxt('data/BETA_GEM.order34.dat').T

w = readmultispec(filename, reform=True, quiet=True)
      # "".join(self.header['WAT2_0*'].values()).split("spec")
      #w = airtovac(gg['wavelen'][orders])
print(w[33])
