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

