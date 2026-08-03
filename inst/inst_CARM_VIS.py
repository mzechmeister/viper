#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file was written by Naïm Teriitehau-Martin. For any issues, please contact : naimteriitehau@gmail.com 
"""

import numpy as np
import os
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.constants import c

from inst.template import read_tpl
from inst.FTS_resample import resample, FTSfits 

'''
If creating or using a template with viper : oset = '0:56'
If using a template from SERVAL : oset = '10:52'
'''
oset = '0:56'    # Relative orders that are analyzed (Here we keep only the 20th to 50th available orders. For CARM_VIS, the absolute orders start at 118 [low wavlengths] and end at 58 [high wavelengths])
iset='400:3600'    # Keep the pixels between X and Y ; Do not analyze the ones outside of that range -- Chosen arbitrarily and can be changed by the user (goes up to 0:4000)
fix = 'wave ip'    # CARMENES is stabilized : need to -fix wave 

# Convert FWHM resolution to sigma
ip_guess = {'s' : c.to_value(u.km/u.s)/(94_600*2*np.sqrt(2*np.log(2)))}  # speed of light [km/s] divided by (spectrograph_resolution x factor) -- Assuming Gaussian IP for factor ; FWHM is defined in terms of speed dv [km/s] (which is ~= delta(ln(wavelen)))  

# Location of CARMENES Spectrograph -- Obtained from the data headers (hdr keys 'HIERARCH CAHA TEL GEOLAT' and 'HIERARCH CAHA TEL GEOLON' and 'HIERARCH CAHA TEL GEOELEV')
location = carmenes = EarthLocation.from_geodetic(lat=37.2236*u.deg, lon=-2.54625*u.deg, height=2168.*u.m)



def Spectrum(filename='', order=None, targ=None):

    hdu = fits.open(filename, ignore_blank=True)
    hdr = hdu[0].header
    
    dateobs = hdr.get('DATE-OBS')   # UTC datetime at observation start (ISOT format)
    exptime_tmean = hdr.get('HIERARCH CARACAL TMEAN') * u.s  # Flux-weighted midpoint of exposure (more accurate than the exptime)
    
    # If target not specified while calling Spectrum() : Define target from data file header info
    ra = hdr.get('RA', np.nan)     # RightAscension of target [degrees]
    dec = hdr.get('DEC', np.nan)    # Declination of target [degrees]
    dec = (dec+90) % 180 - 90   #Declination needs to be between -90 and 90 deg, but some SERVAL tpl headers give values way bigger than that so we fix it
    targdrs = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    if not targ: targ = targdrs
    
    
    # Apply barycentric RV correction halfway through exposure (offers the least error?)
    midtime = Time(dateobs, format='isot', scale='utc') + exptime_tmean
    
    berv = targ.radial_velocity_correction(obstime=midtime, location=carmenes)  #Barycentric Earth RV : This is the RV correction to apply
    berv = berv.to_value(u.km/u.s)
    bjd = midtime.tdb    # Convert midtime scale from utc to Barycentric Dynamical Time
    
    
    # Read file data to obtain spec, wavelen, err
    spec = hdu['SPEC'].data    # Flux for spectrum
    wavelen = hdu['WAVE'].data    # Vacuum wavelength for each pixel [angstrom]
    err_spec = hdu['SIG'].data    # Error estimate for flux
    
    
    # If specific order is selected, only keep that order
    if order is not None:
        wavelen, spec, err_spec = wavelen[order], spec[order], err_spec[order]
    
    # Build pixel array and bad pixel map
    pixel = np.arange(spec.size)    # As many pixels as there are values in spec (Including bad pixels)
    flag_pixel = 1 * np.isnan(spec)  # Bad pixels (spec value undefined) have a value of 1
    
    return pixel, wavelen, spec, err_spec, flag_pixel, bjd, berv



def Tpl(tplname, order=None, targ=None):
    '''
    Tpl should return barycentric corrected wavelengths.
    read_tpl(tplname) reads wavelen and flux from template ; applies necessary corrections depending on type of file (barycentric correc, vacuum correc, ...)
    '''
    wavelen, spec = read_tpl(tplname, inst=os.path.basename(__file__), order=order, targ=targ)

    return wavelen, spec



# If using a cell (default = No cell) -- CARMENES does not use a cell so we should not need this
def FTS(ftsname='None', dv=100):
    '''FTSFits() reads FTS of the cell and obtains wavenumber and flux (f) from data headers
    Converts wavenumbers [cm] to wavelengths (w) [angstrom], also inverts w and f arrays ( [::-1] ) so that wavelengths are in ascending order
    resample() returns w, f, uj= np.arange(ln(w)) with step dv/c, iod_j=flux interpolated from uj
    '''
    return resample(*FTSfits(ftsname), dv=dv)



# For template creation
def write_fits(wave_tpl_all, spec_tpl_all, err_tpl_all, list_files, file_out):
    '''
    When creating a template : All orders seem to give good results except for 56+
    --> Orders 43, 44 scatter a bit more than the rest
    '''
    
    # Get data from header of the first fits file
    file_in = list_files[0]
    hdu = fits.open(file_in, ignore_blank=True)
    
    hdr = hdu[0].header
    
    # grab shape of data arrays
    spec = hdu['SPEC'].data    # Flux for spectrum
    dim = spec.shape
    
    # set everything to nan by default (for orders that are not in the tpl)
    wavelen = np.full(dim, np.nan)   # Vacuum wavelength for each pixel [angstrom]
    spec = np.full(dim, np.nan)      # Flux for spectrum
    err_spec = np.full(dim, np.nan)  # Error estimate for flux
    
    # write template data to the file
    for o in list(spec_tpl_all.keys()):
        wavelen[o] = wave_tpl_all[o]
        spec[o] = spec_tpl_all[o]
        err_spec[o] = err_tpl_all[o]
        
    hdu['SPEC'].data = spec   # Flux for spectrum
    hdu['WAVE'].data = wavelen    # Vacuum wavelength for each pixel [angstrom]
    hdu['SIG'].data = err_spec    # Error estimate for flux
    

    ### Add header cards
    hdr['DATE-TPL'] = (Time.now().isot , 'UT date of template creation')
    
    for n, obsname in enumerate(list_files):
        hdr[f'HIERARCH VIPER COADD FILE {n}'] = os.path.basename(obsname)    # name of coadded files
        
    hdr['HIERARCH VIPER COADD COMIN'] = (min(spec_tpl_all.keys()), 'minimum coadded order')
    hdr['HIERARCH VIPER COADD COMAX'] = (max(spec_tpl_all.keys()), 'maximum coadded order')
    
    hdr['HIERARCH VIPER COADD NUM'] = (len(list_files), 'number of coadded spectra')
    
    
    # hdu.writeto needs to be fixed because some cards are in the wrong order
    ## --> The fix puts EXTNAME(='SIG','WAVE','SPEC') after PCOUNT and GCOUNT instead of before
    hdu.writeto(file_out+'_tpl.fits', output_verify='silentfix', overwrite=True)
    hdu.close()



def get_drift(filepath=''):
    '''
    gather FP drift from the FITS headers ; Requires to re-open every file
    Would be better if this was included in Spectrum() but it does not cost a lot of computing time
    (should take about 30s to run this for 1000 files)
    '''
    
    hdu = fits.open(filepath, ignore_blank=True)
    hdr = hdu[0].header
    
    FP_drift = hdr.get('HIERARCH CARACAL SERVAL FP RV', hdr.get('HIERARCH CARACAL DRIFT FP RV', 0))
    e_drift = hdr.get('HIERARCH CARACAL SERVAL FP E_RV', hdr.get('HIERARCH CARACAL DRIFT FP E_RV', np.inf))

    return FP_drift, e_drift




def driftcorr(tag, obsnames, orders, chunks, tellshift_option=False):

            ''' drift_mes is a 1D array (1 value per file, the same for every order) '''
            drift_mes = np.full(len(obsnames), np.nan)
            e_drift_mes = np.full_like(drift_mes, np.nan)
            
            # initialize rv arrays (rows = files /// columns = orders, chunks)
            rv = np.full((len(obsnames), chunks*len(orders)), np.nan)
            e_rv = np.full_like(rv, np.nan)
            
            # Also use telluric shift for drift correction (e.g. in case drift could not be measured and drift_mes is empty)
            # Currently only for CARMENES
            if tellshift_option:
                # initialize tellshift array (rows = files // columns = orders, chunks)
                tellshift = np.full((len(obsnames), chunks*len(orders)), np.nan)
                e_tellshift = np.full_like(tellshift, np.nan)
             


            # read the .rvo.dat file produced by viper, to gather the data
            rvodata = np.genfromtxt(tag+".rvo.dat", dtype=None, delimiter=" ", names=True, deletechars="~!@#$%^&*()=+~\|]}[{';: /?.>,<", encoding=None)

            # gather data from files
            bjd, RV, e_RV, berv, filename = rvodata['BJD'], rvodata['RV'], rvodata['e_RV'], rvodata['BERV'], rvodata['filename']


            # get the header names for the rv and e_rv columns
            headers = rvodata.dtype.names
            rvcols = [col for col in headers if 'rv' in col and 'e_rv' not in col]    # ['rv0', 'rv1', 'rv2', ...]
            e_rvcols = [col for col in headers if 'e_rv' in col]    # ['e_rv0', 'e_rv1', 'e_rv2', ...]
        
            
            # fill rv array (1 value per chunk) with data from rvo file
            for o in range(len(rvcols)):
                rv[:, o] = rvodata[rvcols[o]]
                e_rv[:, o] = rvodata[e_rvcols[o]]
                
                

            # gather the drift_mes for each file
            ''' 
            With this implementation, drift_mes can only be gathered after the rvo file is written.
                (this means that the RVs returned by viper while it is running, will not be corrected)
            This could be modified by calling get_drift() directly from viper.py (run it 1x per file with the fit_chunk() loop)
                or by returning drift_mes directly from the Spectrum() function
            '''
            for n in range(len(bjd)):    # file number X has index n = (X-1) (file number starts at 1 but index starts at 0)
                drift_mes[n], e_drift_mes[n] = get_drift(obsnames[n])

            
            # If tellshift is ON, mix drift_mes with tellshift
            if tellshift_option:
                
                '''
                Gather the telluric shift :
                Reads the .par.dat file produced by viper, to gather the tellshift and e_tellshift
                '''
                
                # read tellshift from the .par.dat file (after it has been closed)
                pardata = np.genfromtxt(tag+'.par.dat', dtype=None, delimiter=' ', names=True)
                headers = pardata.dtype.names

                # get name of columns for tellshift and e_tellshift
                atmcol = [col for col in headers if 'atm' in col]
                tshiftcol, e_tshiftcol = atmcol[-2], atmcol[-1]    # tellshift and e_tellshift are always last place in par.atm
                print(f'atmcols are {atmcol}')
                
                i_o = pardata['order'] - np.nanmin(orders)    # index of order
                orchunk_index = i_o * chunks + pardata['chunk']    # array-version equivalent of (i_o*chunks+ch) : index of each (order, chunk)
                

                # list of tuples : (tellshift, e_tellshift, n, [i_o*chunks+ch]) ; 1 tuple for each successfully fitted chunk (=1 line from par.dat)
                ziplist = list(zip(pardata[tshiftcol], pardata[e_tshiftcol], pardata['n'], orchunk_index))
                for zipchunk in ziplist:
                    n = zipchunk[2] - 1    # index of file number (file number starts at 1 so we subtract)
                    k = zipchunk[3]    # i_o*chunks + ch
                    
                    # /!\ self.tellshift is expressed in units of [km/s] : Convert it to [m/s]
                    tellshift[n][k] = zipchunk[0] * 1000
                    e_tellshift[n][k] = zipchunk[1] * 1000
                
                
                ''' 
                First, mix drift_mes with weighted mean tellshift (1 drift value per file)
                Then compuite self.tellshift_wmean (1 value per file)
                Then bring self.tellshift_wmean to same median (and same sign) as self.drift_mes (subtract it from RV to get the drift corrected RV)
                Finally, compute the self.drift_mix array (mix of drift_mes and tellshift_wmean)
                '''
                
                ''' Run this 1x after the loop, after computing the RVs (and the tellshifts) for all orders of all files'''
                
                #### remove the median offset from each order of tellshift?
                # This is not currently implemented, but the value of the (weighted mean) tellshift has a different offset for each order, which may need to be removed?
                
                #weighted mean over all orders (1 value per file)
                tellshift_wmean, e_tellshift_wmean = wsem(np.nan_to_num(tellshift, nan=0, posinf=np.inf, neginf=-np.inf), e=np.nan_to_num(e_tellshift, nan=np.inf, posinf=np.inf, neginf=-np.inf), axis=1)
                
                ## if only 1 order could be analyzed, then by default tellshift_wmean is the value of tellshift for that order, and e_tellshift_wmean = 0. 
                # In that case, we replace e_tellshift_wmean with e_tellshift for that order.
                e_tellshift_wmean[np.nansum(np.isfinite(tellshift), axis=1)==1] = np.nansum(e_tellshift, axis=1)[np.nansum(np.isfinite(tellshift), axis=1)==1]

                
                ### tellshift and FP drift have opposite sign : tellshift needs to be ADDED to the RVs while FP drift needs to be SUBTRACTED
                ### Change the sign of tellshift to match the sign of drift_mes : tellshift now needs to be SUBTRACTED from the RVs to apply drift correction
                tellshift_wmean *= -1
                med_tellshift = np.nanmedian(tellshift_wmean)
                med_drift_mes = np.nanmedian(drift_mes)
                
                # bring tellshift and drift_mes to the same median for proper mixing (remove any offset between the two)
                med_drift = med_drift_mes
                tellshift_wmean += med_drift - med_tellshift
                
                
                # drift_mix is a 1D array (1 value per file) because drift_mes is 1D
                drift_mix = np.full_like(drift_mes, np.nan)
                e_drift_mix = np.full_like(drift_mix, np.nan)
                
                
                ### Compute weighted mean drift
                
                # remove nans for proper computation
                e_drift_mes[e_drift_mes==0] = np.inf
                e_drift_mes = np.nan_to_num(e_drift_mes, nan=np.inf, posinf=np.inf, neginf=np.inf)
                
                e_tellshift_wmean[e_tellshift_wmean==0] = np.inf
                e_tellshift_wmean = np.nan_to_num(e_tellshift_wmean, nan=np.inf, posinf=np.inf, neginf=np.inf)
                
                w_drift_mes = 1/(e_drift_mes**2)
                w_tellshift = 1/(e_tellshift_wmean**2)
                w_mix = np.full_like(w_drift_mes, np.nan)    # extra weight that will depend on how reliable drift_mes and tellshift are
                
                # if drift (FP or tshift) not finite then set weight to 0 
                drift_mes_isfinite = np.isfinite(drift_mes)
                tellshift_wmean_isfinite = np.isfinite(tellshift_wmean)
                
                
                for n in range(len(drift_mes)):
                    # If both drift_mes and tellshift are ok, then mix together
                    if drift_mes_isfinite[n] and tellshift_wmean_isfinite[n]:
                        diff_tellshift = np.nan_to_num(np.abs(tellshift_wmean[n] - med_drift), nan=np.inf, posinf=np.inf)
                        diff_drift_mes = np.nan_to_num(np.abs(drift_mes[n] - med_drift), nan=np.inf, posinf=np.inf)
                        
                        # Final weight : w_mix=1 if only use tellshift (=drift_mes is wrong); w_mix=0 means only use drift_mes (tellshift is wrong)
                        #compare how much drift_mes and tellshift_wmean differ from the median drift
                        w_mix[n] = diff_drift_mes / (diff_tellshift + diff_drift_mes)

                        drift_mix[n] = np.nan_to_num((w_mix[n]*(tellshift_wmean[n]*w_tellshift[n]) + (drift_mes[n]*w_drift_mes[n])*(1-w_mix[n])) / (w_drift_mes[n]*(1-w_mix[n]) + w_tellshift[n]*w_mix[n]), nan=0)
                        e_drift_mix[n] = np.sqrt((e_tellshift_wmean[n]*w_tellshift[n]*w_mix[n])**2 + (drift_mes[n]*w_drift_mes[n]*(1-w_mix[n]))**2) 
                        
                    # If both are bad, then keep a nan value
                    elif (not drift_mes_isfinite[n] and not tellshift_wmean_isfinite[n]):
                            pass
                    
                    # If no drift_mes but tellshift is ok : only use tellshift
                    elif not drift_mes_isfinite[n]:
                        drift_mix[n] = tellshift_wmean[n]
                        e_drift_mix[n] = e_tellshift_wmean[n]
                    
                    # If no tellshift but drift_mes is ok : only use drift_mes
                    elif not tellshift_wmean_isfinite[n]:
                        drift_mix[n] = drift_mes[n]
                        e_drift_mix[n] = e_drift_mes[n]

                drift, e_drift = drift_mix, e_drift_mix
                
                
            # If tellshift is OFF, just use drift_mes
            elif not tellshift_option:
                drift, e_drift = drift_mes, e_drift_mes


            ### Apply drift correction to the RVs
            rv -= np.nan_to_num(drift[:, None], nan=0)    # convert drift to a 2D array of shape (N, 1) with None
            e_rv = np.sqrt(e_rv**2 + np.nan_to_num(e_drift[:,None]**2, nan=0))    # propagate errors ; if drift is nan then it is not corrected

            RV -= np.nan_to_num(drift, nan=0)
            # e_RV does not need to be modified (it is a standard dev)


            # Currently this overwrites the .rvo.dat file. Could try writing in a separate '.rvc' file? 
            rvcfile = open(tag+'.rvc.dat', 'w')
            print('BJD RV e_RV BERV', *sum(zip(rvcols, e_rvcols), ()), 'filename', file=rvcfile)
            for n in range(len(bjd)):
                #print(bjd[n], RV[n], e_RV[n], berv[n], drift[n], e_drift[n], *sum(zip(rv[n], e_rv[n]), ()), filename[n], file=rvcfile)
                print(bjd[n], RV[n], e_RV[n], berv[n], *sum(zip(rv[n], e_rv[n]), ()), filename[n], file=rvcfile)
            
            rvcfile.close()

        
