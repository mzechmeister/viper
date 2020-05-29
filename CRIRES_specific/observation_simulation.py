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

def doppler_shift(wvl, flux, v):
     
     # Shifted wavelength axis
     wlprime = wvl * (1.0 + v/299792458)
     
     # Interpolated flux for new wavelength scale
     nflux = pyasl.intep(wlprime, flux, wvl, boundsError=False)
   
     ## Edge handling
     firsts = []
     # Search for first non-NaN value save indices of
     # leading NaN values
     for i in smo.range(len(nflux)):
          if np.isnan(nflux[i]):
               firsts.append(i)
          else:
               firstval = nflux[i]
          break
    # Do the same for trailing NaNs
     lasts = []
     for i in smo.range(len(nflux)-1,0,-1):
          if np.isnan(nflux[i]):
               lasts.append(i)
          else:
               lastval = nflux[i]
          break
    # Use first and last non-NaN value to
    # fill the nflux array
     nflux[firsts] = firstval
     nflux[lasts] = lastval
     return nflux, wlprime

def gaussian(x, amp, cen, sig):
    """1-d gaussian: gaussian(x, amp, cen, sig)"""
    return (amp / (sqrt(2*pi) * sig)) * exp(-(x-cen)**2 / (2*sig**2))

def PSF_convolve(wave,spec,CenGauss_sigma,SatGauss1_amp,SatGauss2_amp):
     
     dxs = wave[1:] - wave[0:-1]
     lx = len(wave)
     nx = (np.arange(lx, dtype=np.int) - sum(divmod(lx, 2)) + 1) * dxs[0]
    

     gauss_a1 = 1
     gauss_s1 = CenGauss_sigma
     gauss_c1 = 0

     gauss_a2 = SatGauss1_amp
     gauss_s2 = 0.03
     gauss_c2 = 0.05

     gauss_a3 = SatGauss2_amp
     gauss_s3 = 0.04
     gauss_c3 = -0.05
    


     mult_gauss = gaussian(nx,gauss_a1,gauss_c1,gauss_s1)+gaussian(nx,gauss_a2,gauss_c2,gauss_s2)+gaussian(nx,gauss_a3,gauss_c3,gauss_s3)
#     mult_gauss = mult_gauss/np.max(mult_gauss)
    

     nf = len(spec)
     spec_1 = np.concatenate((np.ones(nf) * spec[0], spec, np.ones(nf) * spec[-1]))
     result = np.convolve(spec_1, mult_gauss, mode="same")[nf:-nf]

     return result

def model_spec(params, input_data):  
             
     v_star = params['Dop_shift'].value  
     v_gas =  params['Ins_shift'].value
     CenGauss_sigma = params['CenGauss_sigma'].value
     SatGauss1_amp = params['SatGauss1_amp'].value
     SatGauss2_amp = params['SatGauss2_amp'].value
     
     wave = input_data['wave']
     flux_star = input_data['flux_star']
     flux_gas = input_data['flux_gas']
    
    
     Fluxshift_star, Waveshift_star  =doppler_shift(wave, flux_star, v_star)
     Fluxshift_gas, Waveshift_fas  =doppler_shift(wave, flux_gas, v_gas)
  
     product_spec = Fluxshift_star*Fluxshift_gas 
  
    
     convolved_spec = PSF_convolve(wave, product_spec,CenGauss_sigma,SatGauss1_amp,SatGauss2_amp)

     return (convolved_spec)
     
     
def minimize_func(params,input_data, observed_spectrum):
     
     modelled_spec = model_spec(params,input_data)

     return (observed_spectrum-modelled_spec)



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

###################             RVlist READ              ######################

df = pd.read_csv('RV_vals.csv')

RV = df['RV'].values
Dop_vals = df['Dop_shift'].values
ins_vals = df['ins_shift'].values


###############    Wavelength range select (Angstroms)     ####################


wave_select = [13500,13600]

###############           Loop for RV time series        ######################


#for r in range(RV.size):
for r in range(2):
     rv_vals = RV[r]
    
     Dop_shift = np.array([Dop_vals[r]])
     Ins_shift = np.array([ins_vals[r]])

                   
     data_out = pd.DataFrame([])
     observed_spec = pd.DataFrame([])
          
     sigma_gauss1 = 0.01
     
     amp_gauss2 = 0.3
     
     amp_gauss3 = 0.4


     
     
     
     Dop_in = Dop_shift-5
     Ins_in = Ins_shift-5
     sigma1_in = sigma_gauss1 -sigma_gauss1/2
     amp2_in = amp_gauss2 -amp_gauss2/2
     
     amp3_in = amp_gauss3 -amp_gauss3/2

     
     
     
     
     Dop_max = Dop_shift*2
     Ins_max = Ins_shift*2
     sigma1_max = sigma_gauss1*2
     amp2_max = amp_gauss2*2

     amp3_max = amp_gauss3*2
     
     Dop_min = Dop_shift/2
     Ins_min = Ins_shift/2
     sigma1_min = sigma_gauss1/2
     amp2_min = amp_gauss2/2
     
     amp3_min = amp_gauss3/2
     
     
     
     param_in = pd.DataFrame({'Dop_shift':Dop_shift,'Dop_in':Dop_in,'Dop_max':Dop_max,
                           'Dop_min':Dop_min,'Ins_shift':Ins_shift,'Ins_in':Ins_in,
                           'Ins_max':Ins_max,'Ins_min':Ins_min,
                           'sigma1_in':sigma1_in,'sigma1_max':sigma1_max,'sigma1_min':sigma1_min,'sigma_gauss1':sigma_gauss1,
                           'amp2_in':amp2_in,'amp2_max':amp2_max,'amp2_min':amp2_min, 'amp_gauss2':amp_gauss2,
                           'amp3_in':amp3_in,'amp3_max':amp3_max,'amp3_min':amp3_min,'amp_gauss3':amp_gauss3,

                                                                  })
         
###############      Loop for each of the Shifts in RV time series      ######################
         
     for i in range (param_in.Dop_shift.values.size):
          
                         
          wave_split=np.arange(wave_select[0],wave_select[-1],5)
              
          
#          for o in range(wave_split.size-1):
          for o in range(1):
#               o=360
               
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
               

          
               
               w_t = FTS_wavelin
               f_s = Synth_normflux
               f_i = FTS_normflux
     
          
               f_ts, w_ts = doppler_shift(w_t, f_s, param_in['Dop_shift'][i])
               f_is, w_is = doppler_shift(w_t, f_i, param_in['Ins_shift'][i])

          
     #         fig,ax = plt.subplots()
     #         
     #         ax.plot(w_t[0], f_s[0],'b',label='original_Star')
     #         ax.plot(w_ts, f_ts,'r',label='Doppler shifted')
     #         ax.legend()
     #         
     #         fig,ax = plt.subplots()
     #         
     #         ax.plot(w_t[0], f_i[0],'b',label='original_iodine')
     #         ax.plot(w_ts, f_is,'r',label='instrument shifted')
     #         ax.legend()
              
              
               pdt_shift =  f_ts*f_is
                  
                  
               shift_g = PSF_convolve(w_t, pdt_shift, param_in['sigma_gauss1'][i],
                                      param_in['amp_gauss2'][i],param_in['amp_gauss3'][i])
               spec_noise = shift_g + np.random.normal(0.0, 0.01, shift_g.size)


              
               observed_spec = observed_spec.append(pd.DataFrame({'wave':w_t,'spec':spec_noise}))
               
               
               
     rootout = r'/home/sireesha/Desktop/CRIRES/stage4/data/output/observed_spec/'
     observed_spec.to_csv(rootout+'test'+'RV'+str(r)+'.csv')
