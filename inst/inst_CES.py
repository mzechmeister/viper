import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

from .airtovac import airtovac

from .FTS_resample import resample, FTSfits
from pause import pause

location = lasilla = EarthLocation.of_site('lasilla')

pg = {'s': 300_000/200_000/ (2*np.sqrt(2*np.log(2))) }   # convert FHWM resolution to sigma

# see https://github.com/mzechmeister/serval/blob/master/src/inst_FIES.py
#/run/media/zechmeister/5748244b-c6d4-49bb-97d4-1ae0a8ba6aed/data/disk2/zechmeister/CES/reducedv4/ZET1RET/

def Spectrum(filename, o=None, targ=None, chksize=4000):
    with open(filename) as myfile:
        hdr = [next(myfile) for x in range(21)]
    # PX#   WAVELENGTH          FLUX           ERROR         MASK (0/1/6)
    x, w, f, e_f, m = np.genfromtxt(filename, skip_header=21).T
    w = airtovac(w)
    if o is not None:
        o = slice(o*chksize, (o+1)*chksize)
        x, w, f, e_f, m = x[o], w[o], f[o], e_f[o], m[o]

    if 1:
        # stitching
        ax = 0*f + 1    # default pixel size
        ax[576] = 0.963
        ax[4+512*1] = 0.986
        ax[4+512*2-1] = 0.981
        ax[4+512*2] = 0.961
        ax[4+512*3-1] = 0.993
        ax[4+512*4] = 0.987
        ax[4+512*5] = 0.992
        ax[4+512*6-1] = 0.986
        ax[4+512*7] = 0.971
        x = np.cumsum(ax) - 1

    b = 1 * np.isnan(f) # bad pixel map
    #b[f>1.5] |= 2 # large flux
    b[2175:2310] |= 4  # grating ghost on spectrum, CES.2000-08-13T073047.811, stationary?

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
    # sc = SkyCoord(ra=ra, dec=de, unit=(u.hourangle, u.deg), pm_ra_cosdec=pmra*u.mas/u.yr, pm_dec=pmde*u.mas/u.yr)
    midtime = Time(dateobs, format='isot', scale='utc') + exptime * u.s
    berv = targ.radial_velocity_correction(obstime=midtime, location=lasilla)  
    berv = berv.to(u.km/u.s).value  
    bjd = midtime.tdb
    return x, w, f, b, bjd, berv

def Tpl(tplname, o=None, targ=None):
    if tplname.endswith('.dat'):
        # echelle template
        w, f, b, bjd, berv = Spectrum(tplname, targ=targ)
        w *= 1 + berv/3e5
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




