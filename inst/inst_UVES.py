import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.constants import c

from .readmultispec import readmultispec
from .airtovac import airtovac

from .FTS_resample import resample, FTSfits

#location = tls = EarthLocation.from_geodetic(lat=50.980111*u.deg, lon=11.711167*u.deg, height=342*u.m)
location = paranal = EarthLocation.of_site('paranal')

oset = '5:16'
iset = '400:2000'

ip_guess = {'s': 300_000/220_000/ np.sqrt(np.log(256))}   # convert FHWM resolution to sigma

from pause import *
def Spectrum(filename, order=None, targ=None):
    offset = 0
    with open(filename) as hdu:
        while 1:
            line = next(hdu)
            if line.startswith(' 10001'): break
            offset += 1
            if 'ObsDate' in line: date = line.split()[-1]
            #if 'ObsTime' in line: ut = line.split()[-1]  # no seconds!
            if 'ObsTime' in line: ut = float(line.split()[4])
            if 'ExpTime' in line:
                exptime = float(line.split()[-1])

    hh = int(ut)
    mm = int((ut-hh)*60)
    ss = ((ut-hh)*60-mm) * 60
    dateobs = "%s-%s-%sT%s:%s:%s"% (date[-4:], date[2:4], date[:2], hh, mm, ss)

    ra = 0.0
    de = 0.0

    targdrs = SkyCoord(ra=ra, dec=de, unit=(u.deg,u.deg))
    if not targ: targ = targdrs
    midtime = Time(dateobs, format='isot', scale='utc') + exptime/2 * u.s
    berv = targ.radial_velocity_correction(obstime=midtime, location=paranal)
    berv = berv.to(u.km/u.s).value
    bjd = midtime.tdb
 
    d = np.genfromtxt(filename, skip_header=offset+4096*order, dtype='i,f,f,f,i', names='ordpix,wave,flux,e_flux,flag', max_rows=4096)
    w, f, e = d['wave'], d['flux'], d['e_flux']
    w = airtovac(w)

    # fix start guess for some nights (quick hack for the custom reduction)
    if '24052000' in filename: w *= (1-1/3e5)
    if '06082002' in filename: w *= (1+2/3e5)
    if '11072002' in filename: w *= (1+2/3e5)
    if '18092002' in filename: w *= (1+2/3e5)
    if '19092002' in filename: w *= (1+2/3e5)

    x = np.arange(f.size)
    b = 1 * np.isnan(f) # bad pixel map
    #b[f>1.5] |= 2 # large flux
    #b[(5300<w) & (w<5343)] |= 4  # only for HARPS s1d template (this order misses)
    # TLS spectra have a kink in continuum  at about 1700
    # Also the deconv could have a bad wavelength solution.

    return x, w, f, e, b, bjd, berv

def Tpl(tplname, order=None, targ=None):
    '''Tpl should return barycentric corrected wavelengths'''
    if tplname.endswith('.model'):
        # echelle template
        x, w, f, e, b, bjd, berv = Spectrum(tplname, order=order, targ=targ)
        w *= 1 + (berv*u.km/u.s/c).to_value('')   # *model already barycentric corrected (?)
    elif tplname.endswith('_s1d_A.fits'):
        hdu = fits.open(tplname)[0]
        f = hdu.data
        h = hdu.header
        w = h['CRVAL1'] +  h['CDELT1'] * (1. + np.arange(f.size) - h['CRPIX1'])
        w = airtovac(w)
    elif tplname.endswith('_tpl.fits'):
        hdu = fits.open(tplname)[1]
        #from pause import pause; pause()
        w = np.exp(hdu.data['lnwave'])
        f = hdu.data['flux']
    elif tplname.endswith('.dat'):
        import inst.inst_CES
        x, w, f, e, b, bjd, berv = inst.inst_CES.Spectrum(tplname, targ=targ)
        w *= 1 + berv/3e5
    elif tplname.endswith('.ddd'):
        x, w, f, e, b, bjd, berv = Spectrum(tplname, order=order, targ=targ)  # o =17 for CES
        w *= 1 + (berv*u.km/u.s/c).to_value('')
        from scipy.interpolate import CubicSpline
        uj = np.linspace(np.log(w[0]), np.log(w[-1]), 5*w.size)
        fj = CubicSpline(np.log(w), f)(uj)
        w,f = np.exp(uj), fj
    else:
        # long 1d template
        hdu = fits.open(tplname)
        w = hdu[1].data.field('Arg')
        f = hdu[1].data.field('Fun')
        w = airtovac(w)

    return w, f


def FTS(ftsname='lib/UVES/uves_i2_70_wn_cor.fits.gz', dv=100):
    # https://www.eso.org/sci/facilities/paranal/instruments/uves/tools/uves_i2_70_wn_cor.fits.gz
    w,f = FTSfits(ftsname)
    return resample(w, f, dv=dv)

