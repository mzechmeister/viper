#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 16:34:59 2026
This file was written by Naïm Teriitehau-Martin. For any issues, please contact : naimteriitehau@gmail.com 
"""


import numpy as np
from astropy.io import fits
from astropy.time import Time
import astropy.units as u


class correction:
    def __init__(self, instfile, tpl, *args):

        self.instfile = instfile
        self.tpl = tpl
        self.instname = self.instfile.__name__[10:]    # name of the inst file from inst.inst_NAME
        

    def berv(self, bjd, targ):
        
        inst_location = getattr(self.instfile, 'location', None)
        
        berv = targ.radial_velocity_correction(obstime=bjd, location=inst_location)  #Barycentric Earth RV : This is the RV correction to apply
        
        berv = berv.to_value(u.km/u.s)
        
        ## Could also be modified so that the berv correction is done here? and return the corrected wavelen instead
            ### ---> take as input : lnwave_tpl = np.log(wave_tpl)
            ### ==> then in viper.py when defining S_star just call : np.interp(x, correct.berv(np.log(wave_tpl[order]), bjd, targ), spec_tpl[order]) ?
        #lnwave_tpl += np.log(1+ berv/c)
        # return lnwave_tpl
        
        return berv
    
    
    
    
    # For each instrument, define a function get_drift() which returns the drift and error
    # Drift is applied directly to the RVs at the end, after fitting (for error propagation)
    def drift(self, rv, e_rv, filename):

        # if inst file has a function get_drift() that returns drift (e.g. from fits header)
        try:
            drift, e_drift = self.instfile.get_drift(filename=filename)    # this should only require filename and nothing else (order, targ, ...)?
        except:
            drift, e_drift = 0, 0
        
        rv -= drift
        e_rv = np.sqrt(e_rv**2 + e_drift**2)
        return rv, e_rv
        
        
