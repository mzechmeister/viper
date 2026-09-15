#! /usr/bin/env python3
# Licensed under a GPLv3 style license - see LICENSE

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.io import fits


def barycorr(obs, sa, inst):
    # get instrument location from instrument file 
    inst_location = getattr(inst, 'location', None)
    
    # Barycentric Earth motion
    # try inst_loc, otherwise obs.loc
    berv = obs.targ.radial_velocity_correction(obstime=obs.bjd, location=inst_location)        
    berv = (berv-sa).to_value(u.km/u.s)
        
    return berv
    
def proper_motion(midtime, target):
    # apply proper motion and secular acceleration; use information from simbad

    sa = target.sa * (midtime-Time('J2000.0')).to_value('yr') * u.m/u.s
    target.obstime = midtime
    target = target.apply_space_motion(new_obstime=midtime)
        
    return sa, target
        

    
