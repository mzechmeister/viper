#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 11:08:19 2020

@author: Sireesha Chamarthi
"""

import glob
import pandas as pd
import numpy as np
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
from itertools import chain
from astropy.io import fits


#########################  FTS         ########################################

root_FTS = r'/home/sireesha/Desktop/TLS_model/stage1/data/fts/'

Wavenumber_FTS_in, Flux_FTS_in = pyasl.read1dFitsSpec(root_FTS+'TLS_I2_FTS.fits')
Wavenumber_FTS = [h.byteswap().newbyteorder() for h in Wavenumber_FTS_in]
Flux_FTS0 = [h.byteswap().newbyteorder() for h in Flux_FTS_in]

Wavelength_FTS0=1e8/np.array(Wavenumber_FTS) # scale convertion from wavenumber to wavelength (angstrom) 


Wavelength_FTS= np.fliplr([Wavelength_FTS0])[0]
Flux_FTS = np.fliplr([Flux_FTS0])[0]

FTS_dataframe = pd.DataFrame({'wave':Wavelength_FTS,'flux':Flux_FTS})

#########################  Star deconv  #######################################

root_deconv=r'/home/sireesha/Desktop/TLS_model/stage1/data/beta_gem/deconv/'


Stardeconv_files = glob.glob(root_deconv +"deconv*")
Stardeconv_files.sort()

for j in range(1):
     j=3
     Stardeconv_wave,Stardeconv_flux = np.loadtxt(Stardeconv_files[j],delimiter='  ', skiprows=287,unpack=True)

#########################  Star synthetic #######################################

root_synth=r'/home/sireesha/Desktop/TLS_model/stage1/data/beta_gem/synthetic/'

hdu = fits.open(root_synth+'pepsib.20150409.000.sxt.awl.all6')

Starsynth_wave_in = hdu[1].data.field('Arg')
Starsynth_flux_in =hdu[1].data.field('Fun')

Starsynth_wave = [h.byteswap().newbyteorder() for h in Starsynth_wave_in]
Starsynth_flux = [h.byteswap().newbyteorder() for h in Starsynth_flux_in]

Starsynth_dataframe = pd.DataFrame({'wave':Starsynth_wave,'flux':Starsynth_flux})



#########################  Star + iodine#######################################


root_obs = r'/home/sireesha/Desktop/TLS_model/stage1/data/beta_gem/'

StarIodine_files = glob.glob(root_obs +"BETA_GEM.*")
StarIodine_files.sort()

cols = ['wave','flux']

for i in range (1):
     i=17
     StarIodine_data= (pd.read_csv(StarIodine_files[i],delimiter = '  ',names=cols))
     StarIodine_wave=(StarIodine_data['wave'].values)
     StarIodine_flux=(StarIodine_data['flux'].values)
     
     
     w2=[]
     f2=[]
     o2=[]
     
     

     FTS_select = FTS_dataframe[(FTS_dataframe['wave']>StarIodine_wave[0])&(FTS_dataframe['wave']<StarIodine_wave[-1])]
     Starsynth_select = Starsynth_dataframe[(Starsynth_dataframe['wave']>StarIodine_wave[0])&(Starsynth_dataframe['wave']<StarIodine_wave[-1])]
     
     plt.figure(1)
     plt.plot(Stardeconv_wave,(Stardeconv_flux+6)/1.01,'g',label='BETA_GEM_deconvolved')

     plt.plot(Starsynth_select['wave'],Starsynth_select['flux']+4,'r',label='BETA_GEM_synthetic')
#     plt.plot(FTS_select['wave'],FTS_select['flux']+2,'k',label='FTS')
     plt.plot(StarIodine_wave*(9E1+3/3E5),StarIodine_flux/1.04,'b',label='BETA_GEM+Iodine')
     plt.legend()
     
     
     plt.figure(2)
     plt.plot(Stardeconv_wave*(1+(0.1/3E5)),Stardeconv_flux*1.04,'g',label='BETA_GEM_deconvolved')

     plt.plot(Starsynth_select['wave'],Starsynth_select['flux'],'r',label='PEPSI')
     plt.plot(StarIodine_wave*(1-(9/3E5)),StarIodine_flux/1.04,'b',label='BETA_GEM+Iodine')

     plt.legend()
     
     