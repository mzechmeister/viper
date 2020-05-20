
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 14:54:25 2020

@author: Sireesha Chamarthi
"""
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
    

     nf = len(spec)
     spec_1 = np.concatenate((np.ones(nf) * spec[0], spec, np.ones(nf) * spec[-1]))
     result = np.convolve(spec_1, mult_gauss, mode="same")[nf:-nf]
    

     return result

def model_spec(params, input_data):  
             
     Star_shift = params['Dop_shift'].value  
     Gas_shift =  params['Ins_shift'].value
     CenGauss_sigma = params['CenGauss_sigma'].value
     SatGauss1_amp = params['SatGauss1_amp'].value
     SatGauss2_amp = params['SatGauss2_amp'].value
     
     wave = input_data['wave']
     flux_star = input_data['flux_star']
     flux_gas = input_data['flux_gas']
    
    
     Fluxshift_star, Waveshift_star  =doppler_shift(wave, flux_star, Star_shift)
     Fluxshift_gas, Waveshift_fas  =doppler_shift(wave, flux_gas, Gas_shift)
  
     product_spec = Fluxshift_star*Fluxshift_gas 
  
    
     convolved_spec = PSF_convolve(wave, product_spec,CenGauss_sigma,SatGauss1_amp,SatGauss2_amp)

     return (convolved_spec)
     
     
def minimize_func(params,input_data, observed_spectrum):
     
     modelled_spec = model_spec(params,input_data)

     return (observed_spectrum-modelled_spec)




root = r'/home/sireesha/Desktop/CRIRES/stage4/data/'
all_files = glob.glob(root +"*Data.csv")
all_files.sort()

df = pd.read_csv(root +'RV_vals.csv')

RV = df['RV'].values
Dop_vals = df['Dop_shift'].values
ins_vals = df['ins_shift'].values

#for r in range(RV.size):
for r in range(1):
     rv_vals = RV[r]
    
     Dop_shift = np.array([Dop_vals[r]])
     Ins_shift = np.array([ins_vals[r]])

#     for f in range(len(all_files)): 
     for f in range(1):  
     
     
          df = pd.read_csv(all_files[f])
          
          
          df['wave'] = df['wave']*10000
          w_t = np.array_split(df['wave'].values, 450)
          f_s = np.array_split(df['flux_star'].values,450)
          f_i = np.array_split(df['flux_fts'].values, 450)
          
          w_t = np.array(w_t)
          f_s = np.array(f_s)
          f_i = np.array(f_i)
          
          
          
          
          data_out = pd.DataFrame([])
          
          
          
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
         
         
          for i in range (param_in.Dop_shift.values.size):
              
          
          
          
          
          #w_t[0] = w_t[0][0::3]
          #f_s[0] = f_s[0][0::3]
          #f_i[0] = f_i[0][0::3]
          
#              for o in range(w_t.shape[0]):
              for o in range(1):
#                  o=360
          
                  f_ts, w_ts = doppler_shift(w_t[o], f_s[o], param_in['Dop_shift'][i])
                  f_is, w_is = doppler_shift(w_t[o], f_i[o], param_in['Ins_shift'][i])
          
          
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
                  
                  
                  shift_g = PSF_convolve(w_t[o], pdt_shift, param_in['sigma_gauss1'][i],
                                   param_in['amp_gauss2'][i],
                                   param_in['amp_gauss3'][i])
                  
#                  spec_noise = shift_g + np.random.normal(0.0, 0.01, shift_g.size)


              
              
#              fig,ax = plt.subplots()
#              
#              ax.plot(w_t[0], pdt_shift,'b',label='original_Star')
#              ax.plot(w_t[0], shift_g,'r',label='convolved')
#              
#              ax.legend()
#              
              
              
                  
              
                  
              
                  params = Parameters()
                      
                  params.add('Dop_shift', value=param_in['Dop_in'][i],min=param_in['Dop_min'][i],max=param_in['Dop_max'][i])
                  params.add('Ins_shift', value= param_in['Ins_in'][i],min=param_in['Ins_min'][i],max=param_in['Ins_max'][i])
                  params.add('CenGauss_sigma', value= param_in['sigma1_in'][i],min=param_in['sigma1_min'][i],max=param_in['sigma1_max'][i])
                  params.add('SatGauss1_amp', value= param_in['amp2_in'][i],min=param_in['amp2_min'][i],max=param_in['amp2_max'][i])
                  params.add('SatGauss2_amp', value= param_in['amp3_in'][i],min=param_in['amp3_min'][i],max=param_in['amp3_max'][i])
                  
                  
                  
                  
                  
              
                  
                  
                  x= w_t[o]
                  data = shift_g
                  f1=f_s[o]
                  f2=f_i[o]
                  
                  min_func_input = {'wave':x,'flux_star':f1,'flux_gas':f2}
                  
                  result= minimize(minimize_func, params, args= (min_func_input,data),method='least_squares')
                  final = data + result.residual
                  #
                  #
                  #
                  #report_fit(result)
                  #
                  star_pred = result.params['Dop_shift'].value
                  iod_pred = result.params['Ins_shift'].value
                  sig1_pred = result.params['CenGauss_sigma'].value
                  amp2_pred = result.params['SatGauss1_amp'].value
                  amp3_pred = result.params['SatGauss2_amp'].value
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  star_err = result.params['Dop_shift'].stderr
                  iod_err = result.params['Ins_shift'].stderr
                  sig1_err = result.params['CenGauss_sigma'].stderr
                  amp2_err = result.params['SatGauss1_amp'].stderr
                  amp3_err = result.params['SatGauss2_amp'].stderr
                  
                  
                  
                  
                  RV_given = param_in['Dop_shift'][i]-param_in['Ins_shift'][i]
                  RV_pred = star_pred-iod_pred
                  
                  w_start = w_t[o][0]
                  w_end = w_t[o][-1]
                  w_mean = np.mean(w_t[o])
          
          
                  
                  
                  data_out = data_out.append(pd.DataFrame({'Dop_shift':param_in['Dop_shift'][i],'Ins_shift':param_in['Ins_shift'][i],
                                                           'star_pred':star_pred,'iod_pred':iod_pred,'star_err':star_err,'iod_err':iod_err,
                                                           'RV_given':RV_given,'RV_pred':RV_pred,'w_start':w_start,'w_end':w_end,'w_mean':w_mean,
                                                           'sigma_gauss1':sigma_gauss1,'sig1_pred':sig1_pred,'sig1_err':sig1_err,
                                                           'amp_gauss2':amp_gauss2,'amp2_pred':amp2_pred,'amp2_err':amp2_err,
                                                           'amp_gauss3':amp_gauss3,'amp3_pred':amp3_pred,'amp3_err':amp3_err,
     
                                                           
                                                           }, index=[0]), ignore_index=True)
              
            
               
                    
          rootout = '/home/sireesha/Desktop/CRIRES/stage4/data/output'
          
          data_out.to_csv(rootout+'/'+'change_nonoise3'+str(f)+'RV'+str(r)+'.csv')
