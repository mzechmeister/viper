import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

from .airtovac import airtovac

from .FTS_resample import resample, FTSfits
from pause import pause


# see https://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py
#/run/media/zechmeister/5748244b-c6d4-49bb-97d4-1ae0a8ba6aed/data/disk2/zechmeister/CES/reducedv4/ZET1RET/

def Spectrum(filename='data/TLS/other/BETA_GEM.fits', o=None):
    with open(filename) as myfile:
        hdr = [next(myfile) for x in range(21)]
    # PX#   WAVELENGTH          FLUX           ERROR         MASK (0/1/6)
    x, w, f, e_f, m = np.genfromtxt(filename, skip_header=21).T
    w = airtovac(w)

    b = 1 * np.isnan(f) # bad pixel map
    #b[f>1.5] |= 2 # large flux
    #b[(5300<w) & (w<5343)] |= 4  # only for HARPS s1d template (this order misses)

    dateobs = hdr[2].split()[-1]
    exptime = float(hdr[4].split()[-1])
    ra = '03:17:46.1632605674'
    de = '-62:34:31.154247481'
    ra = '03:18:12.8185412558'
    de = '-62:30:22.917300282'
    pmra = 1331.151
    pmde = 648.523
    #from pause import pause; pause()
    #SkyCoord.from_name('M31', frame='icrs')
    lasilla = EarthLocation.of_site('lasilla')
    sc = SkyCoord(ra=ra, dec=de, unit=(u.hourangle, u.deg), pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmde*u.mas/u.yr)
    midtime = Time(dateobs, format='isot', scale='utc') + exptime * u.s
    berv = sc.radial_velocity_correction(obstime=midtime, location=lasilla)  
    berv = berv.to(u.km/u.s).value  
    bjd = midtime.tdb.jd
    return w, f, b, bjd, berv

def Tpl(tplname, o=None):
    if tplname.endswith('.dat'):
        # echelle template
        w, f, *b = Spectrum(tplname)
    elif tplname.endswith('_s1d_A.fits'):
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


def FTS(ftsname='lib/CES/iodine_50_wn.fits', dv=100):
    print('FTS', ftsname)
    # /run/media/zechmeister/5748244b-c6d4-49bb-97d4-1ae0a8ba6aed/data/mybook/Laptop_disk1/CES_iodine/
    return resample(*FTSfits(ftsname), dv=dv)




