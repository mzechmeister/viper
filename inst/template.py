#! /usr/bin/env python3
# Licensed under a GPLv3 style license - see LICENSE

import os
import numpy as np
import importlib
import astropy.units as u
from astropy.constants import c
from astropy.io import fits
from .airtovac import airtovac

def read_tpl(tplname, inst='inst_TLS.py', order=20, targ='None', wmin=3500, wmax=8000):

    successful_read = 0
    
    if tplname.endswith('_s1d_A.fits') or tplname.endswith('.tpl.s1d.fits'):
        print('read HARPS template', tplname)
        
        hdu = fits.open(tplname)[0]
        spec = hdu.data	
        h = hdu.header
        wave = h['CRVAL1'] +  h['CDELT1'] * (1. + np.arange(spec.size) - h['CRPIX1'])
        if tplname.endswith('_s1d_A.fits'):
            wave = airtovac(wave)
        else:
            wave = np.exp(wave)
        successful_read = 1
    
    elif tplname.endswith('.fits') or tplname.endswith('.model'):
       
        try:
          #  print('read '+inst.split('.')[0][5:]+' template', tplname)
            inst = importlib.import_module('inst.'+str(inst)[:-3])
                       
            pixel, wave, spec, err, flag_pixel, bjd, berv = inst.Spectrum(tplname, order=order, targ=targ)
            if not tplname.endswith('_tpl.model') or not tplname.endswith('_tpl.fits'):
                # apply barycentric correction
                wave *= 1 + (berv*u.km/u.s/c).to_value('')
            successful_read = 1
        except:
            pass        
            
        if not successful_read:        
            try:
               # long 1d template
                #print('read 1D template', tplname)
                hdu = fits.open(tplname)
                wave = hdu[1].data.field('Arg')
                spec = hdu[1].data.field('Fun')
                successful_read = 1
            except:
                pass

    elif 'PHOENIX' in str(tplname):
        '''
        Example:
        --------
        >>> read('lte05100-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')

        The wave file must be in the same directory.

        Phoenix spectra might be used to get absolute RVs, Teff, logg and [Fe/H].
        Also a vsini broadening would be needed (see https://github.com/mzechmeister/serval/blob/a348b4ca77e57b0e9e626f8c9fb147f080cc2418/src/serval.py#L228).
        Of course, the precision with depend on the model (mis-)match of the Phoenix spectra.
        '''
        print('read phoenix spectrum', tplname)
        with fits.open(tplname) as hdulist:
           flux = hdulist[0].data
        with fits.open(os.path.join(os.path.dirname(tplname), 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')) as hdulist:
            wave = hdulist[0].data
        imap = (wmin<wave) & (wave<wmax)

        wave = wave[imap]
        spec = flux[imap]
        successful_read = 1
        
    if not successful_read:
     #   print('Error: Template format is not known for the selected instrument.')
        print('\x1b[0;31;40m' +'Error: Template format is not known for the selected instrument.'+ '\x1b[0m')    
        exit()
        
    return wave, spec
        
        
     

