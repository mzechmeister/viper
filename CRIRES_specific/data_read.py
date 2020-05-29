#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 16:28:53 2020

@author: Sireesha Chamarthi
"""
import os
os.chdir('/home/sireesha/Desktop/CRIRES/stage4/data')
import matplotlib.pyplot as plt
import pandas as pd
import peakutils as pk
from peakutils.plot import plot as pplot
import numpy as np
import numpy.polynomial.polynomial as poly
from numpy import exp, pi, sqrt
from lmfit.models import GaussianModel
from itertools import chain
from PyAstronomy import pyasl
from astropy.io import fits
import pandas as pd
import numpy as np
from PyAstronomy import pyasl
from PyAstronomy.pyaC import pyaErrors as PE
import six.moves as smo
from lmfit import minimize, Parameters, Parameter, report_fit
from lmfit import Minimizer, Parameters, report_fit
from lmfit import Model
from uncertainties import unumpy as unp
import glob
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from scipy.ndimage.interpolation import shift
from PyAstronomy import funcFit as fuf
from numpy import exp, loadtxt, pi, sqrt



###################             FTS READ                   ####################
cols = ['wavenum','flux']

FTS_in = pd.read_csv('CRp_SGC2_FTStmpl-HR0p007-WN5000-10000_Hband.dat',delimiter='\t',names=cols)

FTS_wavelength=1e7/(FTS_in.wavenum.values) # scale convertion from wavenumber to wavelength  

FTS_wavearrIN= np.fliplr([FTS_wavelength])[0]
FTS_fluxIN = np.fliplr([FTS_in.flux.values])[0]
FTS_waveIN = FTS_wavearrIN*10

FTS_df = pd.DataFrame({'wave':FTS_waveIN,'flux':FTS_fluxIN})

###################             Spectrum READ              ####################


Synth_datain = fits.open('lte03800-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')

Synth_wavein = fits.open('WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')

Synth_wave0 = [h.byteswap().newbyteorder() for h in Synth_wavein[0].data]
Synth_flux = [h.byteswap().newbyteorder() for h in Synth_datain[0].data]
Synth_wave = np.array(Synth_wave0)*10

Synth_df = pd.DataFrame({'wave':Synth_wave,'flux':Synth_flux})



###############    Wavelength range select (Angstroms)     ####################


wave_select = [13500,13600]

                         
wave_split=np.arange(wave_select[0],wave_select[-1],5)
    

#for o in range(wave_split.size-1):
for o in range(1):
#     o=3
     
     FTS_range = FTS_df.loc[(FTS_df['wave']>wave_split[o]) & (FTS_df['wave']<wave_split[o+1])]


     FTS_wavearr =FTS_range['wave'].values
     FTS_fluxarr = FTS_range['flux'].values
    

     
     
     Synth_range = Synth_df.loc[(Synth_df['wave']>wave_split[o]) & (Synth_df['wave']<wave_split[o+1])]
     
     Synth_range= Synth_range.sort_values(by=['wave'])
     
     
     Synth_wavearr =Synth_range['wave'].values
     Synth_fluxarr = Synth_range['flux'].values
         
     
     Synthflux_interpolated = pyasl.intep(Synth_wavearr,Synth_fluxarr,FTS_wavearr,boundsError=False)

         
     ## rearraging the wavelength on a linear scale
   
     FTS_wavelin = np.linspace(FTS_wavearr[0],FTS_wavearr[-1],FTS_wavearr.size)
    
     ## fitting a 2nd order polynomial
         
     FTS_polycoefs = poly.polyfit(FTS_wavearr, FTS_fluxarr, 2)
     FTS_polyfit = poly.polyval(FTS_wavelin, FTS_polycoefs)
     
         
     FTS_contflux = FTS_fluxarr/FTS_polyfit    
         
     FTS_normflux = FTS_contflux/np.max(FTS_contflux) ## Normalized flux
     
     
#               plt.figure(1)
#               plt.plot(FTS_wavearr,FTS_fluxarr)
#               plt.plot(FTS_wavelin, FTS_normflux,'b')
     
         
     ## fitting a 2nd order polynomial
     Synth_polycoefs = poly.polyfit(FTS_wavearr, Synthflux_interpolated, 2)
     Synth_polyfit = poly.polyval(FTS_wavelin, Synth_polycoefs)
     
     Synth_contflux = Synthflux_interpolated/Synth_polyfit    
         
     Synth_normflux = Synth_contflux/np.max(Synth_contflux) ## Normalized flux
     
#               plt.figure(2)
#               #plt.plot(FTS_wavearr,Synth_fluxarr)
#               plt.plot(FTS_wavelin, Synth_normflux,'b')
     
     
     ### log scale
     
     wave = np.log(FTS_wavelin)
     