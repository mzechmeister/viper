#! /usr/bin/env python3
## Licensed under a GPLv3 style license - see LICENSE
## viper - Velocity and IP Estimator
## Copyright (C) Mathias Zechmeister and Jana Koehler


import argparse
import glob
import importlib
import os
import time
from collections import defaultdict
import configparser

import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
from astropy.io import fits
import astropy.units as u

from utils.gplot import *
gplot.colors('classic')
gplot2 = Gplot()
from utils.param import Params
from utils.pause import pause

from utils.model import model, model_bnd, IPs, show_model, pade
from utils.targ import Targ
import utils.convert_output as convert_output
try:
    import viper.vpr as vpr
except: 
    import vpr


viperdir = os.path.dirname(os.path.realpath(__file__)) + os.sep

c = 299792.458   # [km/s] speed of light

targ = None
modset = {}   # model setting parameters
insts = [os.path.basename(i)[5:-3] for i in glob.glob(viperdir+'inst/inst_*.py')]

class nameddict(dict):
    """
    Examples
    --------
    >>> nameddict({'a':1, 'b':2})
    {'a': 1, 'b': 2}
    >>> x = nameddict(a=1, b=2)
    >>> x.a
    1
    >>> x['a']
    1
    >>> x.translate(3)
    ['a', 'b']
    """
    __getattr__ = dict.__getitem__

    def translate(self, x):
        return [name for name,f in self.items() if (f & x) or f==x==0]

# flag_pixel map flagging
flag = nameddict(
    ok=       0, # good pixel
    nan=      1, # nan flux in pixel of spectrum
    neg=      2, # significant negative flux in pixel of spectrum (f < -3*f_err < 0 && )
    sat=      4, # too high flux (saturated)
    atm=      8, # telluric at wavelength of spectrum
    sky=     16, # sky emission at wavelength of spectrum
    out=     32, # region outside the template
    clip=    64, # clipped value
    lowQ=   128, # low Q in stellar spectrum (mask continuum)
    badT=   256, # bad corresponding region in the template spectrum
    chunk=  512, # chunk cutting
)


def arg2slice(arg):
    """Convert string argument to a slice."""
    # We want four cases for indexing: None, int, list of ints, slices.
    # Use [] as default, so 'in' can be used.
    if isinstance(arg, str):
        arg = eval('np.s_['+arg+']')
    return [arg] if isinstance(arg, int) else arg

def arg2range(arg):
    return  eval('np.r_['+arg+']')


def SSRstat(vgrid, SSR, dk=1, plot='maybe', N=None):
    '''
    Analyse chi2 peak.

    Parameters
    ----------
    N: Number of data points in the fit. Needed to estimate 1 sigma uncertainty and when SSR are not chi2 values.

    '''
    k = np.argmin(SSR[dk:-dk]) + dk   # best point (exclude borders)
    vpeak = vgrid[k-dk:k+dk+1]
    SSRpeak = SSR[k-dk:k+dk+1] - SSR[k]
    v_step = vgrid[1] - vgrid[0]
    # interpolating parabola a0+a1*x+a2*x**2 (direct solution) through the three pixels in the minimum
    a = np.array([0, (SSR[k+dk]-SSR[k-dk])/(2*v_step), (SSR[k+dk]-2*SSR[k]+SSR[k-dk])/(2*v_step**2)])  # interpolating parabola for even grid
    v = (SSR[k+dk]-SSR[k-dk]) / (SSR[k+dk]-2*SSR[k]+SSR[k-dk]) * 0.5 * v_step

    v = vgrid[k] - a[1]/2./a[2]   # position of parabola minimum
    e_v = np.nan
    if -1 in SSR:
        print('opti warning: bad ccf.')
    elif a[2] <= 0:
        print('opti warning: a[2]=%f<=0.' % a[2])
    elif not vgrid[0] <= v <= vgrid[-1]:
        print('opti warning: v not in [va,vb].')
    else:
        e_v = 1. / a[2]**0.5
        if N:
            # Rescale the variances (parabola) such that the minimum has value N, implying chi^2_red == 1.
            # Then derive one sigma uncertainty from Delta chi^2 = 1.
            SSRmin = SSR[k] + a[0] - a[1] * a[1]/2./a[2] + a[2] * (a[1]/2./a[2]) **2
            e_v *= (SSRmin/N) **0.5


    if (plot==1 and np.isnan(e_v)) or plot==2:
        gplot2.yrange('[*:%f]' % np.max(SSR))
        gplot2(vgrid, SSR-SSR[k], " w lp, vk="+str(vgrid[k])+", %f+(x-vk)*%f+(x-vk)**2*%f," % tuple(a), [v,v], [0,SSR[1]], 'w l t "%f km/s"'%v)
        gplot2+(vpeak, SSRpeak, ' lt 1 pt 6; set yrange [*:*]')
        pause(v)
    return v, e_v, a

    
if __name__ == "__main__" or __name__ == "viper.viper":
    # Print defaults, but do not wrap lines
    argparse.ArgumentDefaultsHelpFormatter._split_lines = lambda self, text, width: text.splitlines()

    # check first the instrument with preparsing
    preparser = argparse.ArgumentParser(add_help=False)
    preparser.add_argument('args', nargs='*')
    preparser.add_argument('-inst', help='Instrument.', default='TLS', choices=insts)
    preparser.add_argument('-config_file', help='Config file and optional section  [None DEFAULT].', nargs='*', type=str)
    preargs = preparser.parse_known_args()[0]
    
    Inst = importlib.import_module('inst.inst_'+preargs.inst)
    FTS = Inst.FTS
    Tpl = Inst.Tpl
    Spectrum = Inst.Spectrum
    Tell = getattr(Inst, 'Tell', None)
    iset = getattr(Inst, 'iset', slice(None))
    oset = getattr(Inst, 'oset')

    # read in default values from config_viper.ini
    configs_inst, configs_user = {}, {}
    config_default = configparser.ConfigParser()
    config_default.read(viperdir+'config_viper.ini')	# default file
    configs_def = dict(config_default['DEFAULT'])		# default values
    if preargs.inst in config_default.sections():
        configs_inst = dict(config_default[preargs.inst])	# instrument values

    if preargs.config_file:
        if len(preargs.config_file) == 1 and not preargs.config_file[0].endswith('.ini'): 
            # use selected section from default .ini file
            if preargs.config_file[0] in config_default.sections():
                configs_user = dict(config_default[preargs.config_file[0]])
        elif len(preargs.config_file) == 2:    
            # Read in values from user #.ini   
            config = configparser.ConfigParser()         
            config.read(preargs.config_file[0])
            if preargs.config_file[1] in config.sections():
                configs_user = dict(config[preargs.config_file[1]])
            else:
                print('WARNING: Declared section is not found in %s. Use DEFAULT values instead.' % preargs.config_file[0])

    parser = argparse.ArgumentParser(description='VIPER - velocity and IP Estimator', add_help=False, formatter_class=argparse. ArgumentDefaultsHelpFormatter)
    argopt = parser.add_argument   # function short cut
    argopt('obspath', help='Filename of observation.', default='data/TLS/betgem/BETA_GEM.fits', type=str)
    argopt('tplname', help='Filename of template.', nargs='?', type=str)
    argopt('-inst', help='Instrument.', default='TLS', choices=insts)
    argopt('-fts', help='Filename of FTS Cell.', default=viperdir + FTS.__defaults__[0], dest='ftsname', type=str)
    argopt('-ip', help='IP model (g: Gaussian, ag: asymmetric (skewed) Gaussian, sg: super Gaussian, bg: biGaussian, mg: multiple Gaussians, mcg: multiple central Gaussians, bnd: bandmatrix).', default='g', choices=[*IPs], type=str)
    argopt('-chunks', nargs='?', help='Divide one order into a number of chunks.', default=1, type=int)
    argopt('-config_file', nargs='*', help='Config file and optional section  [None DEFAULT].', type=str)
    argopt('-createtpl', nargs='?', help='Removal of telluric features (or cell lines) and combination of several observations.', default=False, const=True, type=int)
    argopt('-deg_bkg', nargs='?', help='Number of additional parameters.', default=0, const=1, type=int)
    argopt('-deg_norm', nargs='?', help='Polynomial degree for flux normalisation.', default=3, type=int)
    argopt('-deg_norm_rat', nargs='?', help='Rational polynomial degree of denominator for flux normalisation.', type=int)
    argopt('-deg_wave', nargs='?', help='Polynomial degree for wavelength scale l(x).', default=3, type=int)
    argopt('-demo', nargs='?', help='Demo plots. Use -8 to skip plots 1,2,4).', default=0, const=-1, type=int)
    argopt('-flagfile', help='Use just good region as defined in flag file.', default='', type=str)
    argopt('-infoprec', help='Prints and plots information about precision estimates for the star and the iodine.', action='store_true')
    argopt('-iphs', nargs='?', help='Half size of the IP.', default=50, type=int)
    argopt('-ipB', nargs='*', help='Factor of IP width varation.', type=float, default=[])
    argopt('-iset', help='Pixel range.', default=iset, type=arg2slice)
    argopt('-kapsig', nargs='*', help='Kappa sigma values for the clipping stages. Zero does not clip.', default=[0], type=float)
    argopt('-kapsig_ctpl', help='Kappa sigma values for the clipping of outliers in template creation.', default=0.6, type=float)
    argopt('-look', nargs='?', help='See final fit of chunk with pause.', default=[], const=':100', type=arg2range)
    argopt('-lookfast', nargs='?', help='See final fit of chunk without pause.', default=[], const=':100', type=arg2range)
    argopt('-lookguess', nargs='?', help='Show initial model.', default=[], const=':100', type=arg2range)
    argopt('-lookpar', nargs='?', help='See parameter of chunk.', default=[], const=':100', type=arg2range)
    argopt('-lookres', nargs='?', help='Analyse the residuals.', default=[], const=':100', type=arg2range)
    argopt('-lookctpl', nargs='?', help='Show created template.', default=[], const=':100', type=arg2range)
    #argopt('-nexcl', help='Pattern ignore', default=[], type=arg2range)
    argopt('-molec', nargs='*', help='Molecular specifies; all: Automatic selection of all present molecules.', default=['all'], type=str)
    argopt('-nexcl', nargs='*', help='Ignore spectra with string pattern.', default=[], type=str)
    argopt('-nocell', help='Do the calibration without using the FTS.', action='store_true')
    argopt('-nset', help='Index for spectrum.', default=':', type=arg2slice)
    argopt('-oset', help='Index for order.', default=oset, type=arg2slice)
    argopt('-output_format', nargs='*', help='Format of output files for rvo and par data (dat, fits, cpl).', default=['dat'], dest='oformat', type=str)
    argopt('-oversampling', help='Oversampling factor for the template data.', default=None, type=int)
    argopt('-rv_guess', help='RV guess.', default=1., type=float)   # slightly offsetted
    argopt('-stepRV', help='Step through fixed RVs to find the minimum in the rms (a: (auto) picks the fixed RVs automatically to get close to the minimum; m: (manual) uses fixed range and steps around vguess).', choices=['a', 'm'], type=str)
    argopt('-tag', help='Output tag for filename.', default='tmp', type=str)
    argopt('-targ', help='Target name requested in simbad for coordinates, proper motion, parallax and absolute RV.', dest='targname')
    argopt('-tellshift', nargs='?', help='Variable telluric wavelength shift (one value for all selected molecules).', default=False, const=True, type=int)
    argopt('-telluric', help='Treating tellurics (mask: mask tellurics; sig: downweight tellurics; add: telluric forward modelling with one coeff for each molecule; add2: telluric forward modelling with combined coeff for non-water molecules).', default='', type=str)
    argopt('-tsig', help='(Relative) sigma value for weighting tellurics.', default=1, type=float)
    argopt('-vcut', help='Trim the observation to a range valid for the model [km/s]', default=100, type=float)
    argopt('-wgt', nargs='?', help='Weighted least square fit (employ data error).', default=False, const=True, type=int)
    argopt('-?', '-h', '-help', '--help', help='Show this help message and exit.', action='help')

    parser.set_defaults(**configs_def)
    parser.set_defaults(**configs_inst)
    parser.set_defaults(**configs_user)
   
    parser.set_defaults(kapsig = [float(i) for i in (argopt('--kapsig').default.split(' '))])
    
    args = parser.parse_args()
    globals().update(vars(args))


def fit_chunk(order, chunk, obsname, targ=None, tpltarg=None):
    ####  observation  ####
    pixel, wave_obs, spec_obs, err_obs, flag_obs, bjd, berv = Spectrum(obsname, order=order, targ=targ)

    if telluric == 'mask':
        flag_obs[mskatm(wave_obs) > 0.1] |= flag.atm
    flag_obs[np.isnan(spec_obs)] |= flag.nan

    # select common wavelength range
    lmin = max(wave_obs[iset][0], wave_tpl[order][0], wave_cell[0])
    lmax = min(wave_obs[iset][-1], wave_tpl[order][-1], wave_cell[-1])

    # trim the observation to a range valid for the model
    #  vcut = 100   # [km/s]
    flag_obs[np.log(wave_obs) < np.log(lmin)+vcut/c] |= flag.out
    flag_obs[np.log(wave_obs) > np.log(lmax)-vcut/c] |= flag.out

    # using the supersampled log(wavelength) space with knot index j
    sj = slice(*np.searchsorted(lnwave_j_full, np.log([lmin, lmax])))
    lnwave_j = lnwave_j_full[sj]
    spec_cell_j = spec_cell_j_full[sj]

    ibeg, iend = np.where(flag_obs&1==0)[0][[0, -1]]   # the first and last pixel that is not trimmed
    len_ch = int((iend-ibeg)/chunks)
    ibeg = ibeg + chunk*len_ch
    iend = ibeg + len_ch
    if chunks > 1:
        # divide dataset into chunks
        flag_obs[:ibeg] |= flag.chunk
        flag_obs[iend:] |= flag.chunk

    if flagfile:
        # clip selected pixel ranges in order          
        for msk_i in msk_o[msk_o.order==order]:
            flag_obs[int(msk_i.start):int(msk_i.end)] |= flag.clip        
        # clip selected wavelength regions
        if len(msk_l):   
            msk_wave = lambda x: np.interp(x, msk_w, msk_f)
            flag_obs[msk_wave(wave_obs) > 0.1] |= flag.clip

    if 1:
        # preclip upper outlier (cosmics)
        kap = 6
        p17, smod, p83 = np.percentile(spec_obs[flag_obs==0], [17, 50, 83])
        sig = (p83 - p17) / 2
        flag_obs[spec_obs > smod+kap*sig] |= flag.clip
        # gplot(spec_obs, f', {p17}, {smod}, {p83}, {smod+ kap*sig}')

    # select good pixel
    i_ok = np.where(flag_obs==0)[0]
    pixel_ok = pixel[i_ok]
    wave_obs_ok = wave_obs[i_ok]
    spec_obs_ok = spec_obs[i_ok]

    modset['xcen'] = xcen = np.nanmean(pixel_ok) + 18   # slight offset, then it converges for CES+TauCet
    modset['IP_hs'] = iphs

    if deg_norm_rat:
        # rational polynomial
        modset['func_norm'] = lambda x, par_norm: pade(x, par_norm[:deg_norm+1], par_norm[deg_norm+1:])

    specs_molec = []
    par_atm = parfix_atm = []
    if 'add' in telluric:
        # select present molecules for telluric forward modeling
        specs_molec = np.zeros((0, len(lnwave_j)))
        for mol in specs_molec_all.keys():
            s_mol = slice(*np.searchsorted(wave_atm_all[mol], [lmin, lmax]))
            # bring it to same log(wavelength) scale as cell        
            if specs_molec_all[mol][s_mol] != []:
                spec_mol = np.interp(lnwave_j, np.log(wave_atm_all[mol][s_mol]), specs_molec_all[mol][s_mol])
                specs_molec = np.r_[specs_molec, [spec_mol]]
                # chose just present molecules in wavelength range
                if np.nanstd(spec_mol) > 0.0001:
                    par_atm.append((1, np.inf))
                else:
                   # fix parameter and set it to nan if molecule is not present in order
                    par_atm.append((np.nan, 0))	# fix parameter
            else:
                # set default spectrum if molecule is not present in wavelength range
                specs_molec = np.r_[specs_molec, [lnwave_j*0+1]]
                par_atm.append((np.nan, 0))	# fix parameter 

        if telluric == 'add2' and len(molec) > 1:
            # use combined coeff for all non-water tellurics instead of one for each molecule
            # water tellurics grow with airmass and pwv
            # non-water telluics grow with airmass and depend on seasonal changes
            par_atm = np.asarray(par_atm)
            is_H2O = np.asarray(molec) == 'H2O'

            if any(is_H2O):
                specs_molec = [specs_molec[is_H2O][0], np.nanprod(specs_molec[~is_H2O]*(par_atm[~is_H2O][:, 0]).reshape(-1, 1), axis=0)]
                par_atm = [(1, np.inf), (1, np.inf)]
            else:
                specs_molec = np.nanprod(specs_molec[~is_H2O]*(par_atm[~is_H2O][:, 0]).reshape(-1, 1), axis=0)
                par_atm = [(1, np.inf)]
               
        # add parameter for telluric position shift if selected
        if tellshift:
            par_atm.append((1, np.inf))

    if demo & 1:
        # pre-look raw input
        s_cell = slice(*np.searchsorted(wave_cell, [lmin, lmax]))
        s_tpl = slice(*np.searchsorted(wave_tpl[order], [lmin, lmax]))

        # plot data, template, and iodine with some scaling
        gplot.xlabel('"Vacuum wavelength [Å]"')
        gplot.ylabel('"flux"')
        gplot(wave_cell[s_cell], spec_cell[s_cell]/np.nanmedian(spec_cell[s_cell]), 'w l lc 9 t "cell",', wave_tpl[order][s_tpl], spec_tpl[order][s_tpl]/np.nanmedian(spec_tpl[order][s_tpl]), 'w l lc 3 t "tpl",', wave_obs, spec_obs/np.nanmedian(spec_obs), 'w lp lc 1 pt 7 ps 0.5 t "obs"')
        pause('demo 1: raw input')


    # convert discrete template into a function
    if tplname:
        S_star = lambda x: np.interp(x, np.log(wave_tpl[order]) - np.log(1+berv/c), spec_tpl[order])  # Apply barycentric motion
    else:
        S_star = lambda x: 0*x + 1

    IP = IPs[ip]

    # setup the model
    S_mod = model(S_star, lnwave_j, spec_cell_j, specs_molec, IP, **modset)

    if demo & 2:
        # plot the IP
        gplot.xlabel('"[km/s]"')
        gplot.ylabel('"contribution"')
        gplot(S_mod.vk, S_mod.IP(S_mod.vk), 't "IP model"')
        pause('demo 2: default IP')

    if demo & 4:
       # plot again, now the stellar template can be interpolated
       gplot.xlabel('"Vacuum wavelength [Å]"')
       gplot.ylabel('"flux"')
       gplot(np.exp(lnwave_j), spec_cell_j, S_star(lnwave_j)/np.nanmedian(S_star(lnwave_j)), 'w l lc 9 t "cell", "" us 1:3 w l lc 3 t "tpl"')
       pause('demo 4: stellar template evaluated at lnwave_j')


    # an initial parameter set
    par = Params()

    # a good guess for the stellar RV is needed
  #  par.rv = rv_guess if (tplname or createtpl) else (0, 0)   # else: do not fit for RV
    par.rv = rv_guess if tplname else (0, 0)   # else: do not fit for RV

    # guess for normalization
    norm_guess = np.nanmean(spec_obs_ok) / np.nanmean(S_star(np.log(wave_obs_ok))) / np.nanmean(spec_cell_j)
    par.norm = [norm_guess] + [0]*deg_norm

    if deg_norm_rat:
        # rational polynom
        par.norm += [5e-7] * deg_norm_rat   # a tiny scale hint (zero didn't iterate)
        #par.norm += [ 5e-7**(i+1) for i in range(deg_norm_rat)]   # a tiny scale hint (zero didn't iterate)

    # guess wavelength solution
    par.wave = np.polyfit(pixel_ok-xcen, wave_obs_ok, deg_wave)[::-1]
    parguess = Params(par)

    # guess IP - read in from instrument file
    par.ip = [Inst.ip_guess['s']]
    par.atm = par_atm

    # guess additional background
    if deg_bkg:
        par.bkg = [0] #* deg_bkg

    if demo:
        # disturb guess
        par.norm = parguess.norm = [norm_guess*1.3] + [0]*deg_norm
        # b = par_wave_guess = [wave_ob[0], (wave_obs[-1]-wave_obs[0])/wave_obs.size] # [6128.8833940969, 0.05453566108124]
        par.wave = parguess.wave = [*np.polyfit(pixel[[400, -300]]-xcen-10, wave_obs[[400, -300]], 1)[::-1]] + [0]*(deg_wave-1)
        par.ip = [par.ip[0]*1.5]

    if ip in Inst.ip_guess:
        par.ip = Inst.ip_guess[ip]
    elif ip in ('sg', 'mg', 'asg'):
        par.ip += [2.]   # exponent of super Gaussian
    elif ip in ('ag', 'agr', 'asg'):
        par.ip += [1.]   # skewness parameter (offset to get iterations)
    elif ip in ('bg',):
        par.ip += [par.ip[-1]]   # symmetric biGaussian
    parguess.ip = par.ip

    # set weighting parameter for tellurics
    sig = 1 * err_obs if wgt else np.ones_like(spec_obs)
    if telluric in ('sig', 'add', 'add2'):
        sig[mskatm(wave_obs) < 0.1] = tsig

    if demo & 8:
        # a simple call to the forward model
        # Si_mod = S_mod(pixel_ok, par_rv=0, a=a, b=b, s=s)
        # show the start guess
        S_mod.show(par, pixel_ok, spec_obs_ok, res=False, dx=0.1)
        pause('demo 8: Smod simple call')

    fixed = lambda x: [(pk, 0) for pk in x]
    if demo & 16:
        # A wrapper to fit the continuum
        par_d16 = Params(rv=(rv_guess, 0), norm=[norm_guess], wave=fixed(parguess.wave), ip=fixed(parguess.ip), atm=fixed(par.atm))
        p_norm, _ = S_mod.fit(pixel_ok, spec_obs_ok, par_d16, res=False, dx=0.1, sig=sig[i_ok])
        parguess.norm[0] = p_norm.norm[0]
        pause('demo 16: S_par_norm')

    if demo & 32:
        # A wrapper to fit the wavelength solution
        par_d32 = Params(rv=(rv_guess, 0), norm=fixed(parguess.norm), wave=parguess.wave[:-1]+[1e-15], ip=fixed(parguess.ip), atm=fixed(par.atm), bkg=[(0, 0)])
        p_wave, _ = S_mod.fit(pixel_ok, spec_obs_ok, par_d32, res=False, dx=0.1, sig=sig[i_ok])
        par.wave = p_wave.wave
        pause('demo 32: S_par_wave')

    if demo & 64:
        # fit par_rv, a0 and b simultaneously
        par_d64 = Params(rv=(rv_guess, 0), norm=parguess.norm, wave=par.wave, ip=fixed(parguess.ip), atm=fixed(par.atm), bkg=[(0, 0)])
        params, _ = S_mod.fit(pixel_ok, spec_obs_ok, par_d64, res=False, dx=0.1, sig=sig[i_ok])
        par = Params(params)
        # remove uncertainties
        par = par + dict([(k, v.value) for k,v in par.flat().items()])
        pause('demo 64: S_par_norm_wave_rv')


    if ip in ('sg', 'ag', 'agr', 'bg', 'bnd'):
        # prefit with Gaussian IP
        S_modg = model(S_star, lnwave_j, spec_cell_j, specs_molec, IPs['g'], **modset)

        par1 = Params(par, ip=par.ip[0:1])   # fit only sigma
        par2, _ = S_modg.fit(pixel_ok, spec_obs_ok, par1, sig=sig[i_ok])

        par = par + par2.flat()   # update, but replace first ip par
    par3 = par

    if order in lookguess:
        if demo:
            par_wave_guess = par_wave
            par_norm = [norm_guess]
        params_guess = Params(rv=par.rv, norm=par.norm, wave=parguess.wave, ip=par.ip, atm=parfix_atm, bkg=par.bkg)
        prms = S_mod.show(params_guess, pixel_ok, spec_obs_ok, res=True, dx=0.1)
        pause('lookguess')


    if kapsig[0]:
        # first kappa sigma clipping of outliers
        smod = S_mod(pixel, **par3)
        resid = spec_obs - smod
        resid[flag_obs != 0] = np.nan

        flag_obs[abs(resid) >= (kapsig[0]*np.nanstd(resid))] |= flag.clip
        i_ok = np.where(flag_obs == 0)[0]
        pixel_ok = pixel[i_ok]
        wave_obs_ok = wave_obs[i_ok]
        spec_obs_ok = spec_obs[i_ok]

    if IP == 'bnd':
        # Non parametric fit with band matrix
        # We step through velocity in 100 m/s step. At each step there is linear least square
        # fit for the 2D IP using band matrix.
        S_mod = model_bnd(S_star, lnwave_j, spec_cell_j, params[2], **modset)
        opt = {'x': pixel_ok, 'sig_k': par_ip[0]/1.5/c}
        rr = S_mod.fit(spec_obs_ok, 0.1, **opt)
        fx = S_mod(pixel_ok, 0.1, rr[0])
        ipxj = S_mod.IPxj(rr[0])
        if demo & 2:
            gplot(ipxj, 'matrix w image')

        e_v = np.nan
        if tplname:
            vv = np.arange(-1, 1, 0.1)
            RR = []
            aa = []
            for v in vv:
                rr = S_mod.fit(spec_obs_ok, v, **opt)
                RR.append(*rr[1])
                aa.append(rr[0])
                if 1:
                    print(v, rr[1])

            par_rv, e_v, a = SSRstat(vv, RR, plot=1, N=spec_obs_ok.size)

        best = S_mod.fit(spec_obs_ok, par_rv, **opt)
        fx = S_mod(pixel_ok, par_rv, best[0])
        pause()
        S_mod.show([par_rv, best[0]], pixel_ok, spec_obs_ok, x2=pixel_ok)
        res = spec_obs_ok - fx
        np.savetxt('res.dat', list(zip(pixel_ok, res)), fmt="%s")
        prms = np.nanstd(res) / fx.nanmean() * 100
        if order in look:
            pause()
        return par_rv*1000, e_v*1000, bjd.jd, berv, best[0], np.diag(np.nan*best[0]), prms

    if stepRV in ['a', 'm']:
        # calculating the best RV by going through different fixed RVs
        # finding best value for minimum in rms
        # still some improvement/testing
        rv_range = 0.5
        rms_start = [0, 0, 0]
        rms_all = []
        v_all = []
        rounds = 0       # to make sure, it will not end in an endless loop

        if stepRV == 'm':
            # fix step size in given range
            # otherwise search minimum
            v_grid = np.arange(rv_guess-1, rv_guess+1, 0.1)
            # v_grid = np.arange(-0.1, 0.3, 0.01)
        elif stepRV == 'a':
            v_grid = [rv_guess-rv_range, rv_guess, rv_guess+rv_range]

        while rounds < 20:
            for vv,vguess in enumerate(v_grid):
                if vguess not in v_all:
                    params, e_params = S_mod.fit(pixel_ok, spec_obs_ok, None, par_norm, par_wave_guess, par_ip, par_atm, parfix_rv = vguess, par_bkg=par_bkg, parfix_bkg=parfix_bkg, dx=0.1, sig=sig[i_ok])
                    rms_start1 = np.nanstd(spec_obs_ok - S_mod(pixel_ok, *params))
                    if stepRV == 'a':
                        rms_start[vv] = rms_start1
                    # print('p:', rms_start[vv], p)
                    v_all.append(vguess)
                    rms_all.append(rms_start1)
                else:
                    rms_start[vv] = rms_all[(np.argwhere(np.asarray(v_all)==vguess))[0][0]]

            if (((3 <= np.argmin(rms_all) <= rounds-3) or (abs(v_grid[0]-v_grid[1]) < 0.03)) and (rounds >= 6)) or (stepRV == 'm'):
                # fitting process is done
                rounds = 20
            else:
                # find the position of the current minimum and search further in this direction
                ind = np.argsort(rms_start)
                if ind[0] == 0:
                    v_grid = [v_grid[0]-rv_range, v_grid[0], v_grid[1]]
                elif ind[0] == 1:
                    v_grid = [v_grid[ind[0]], (v_grid[ind[1]]+v_grid[ind[0]])/2, v_grid[ind[1]]]
                elif ind[0] == 2:
                    v_grid = v_grid = [v_grid[1], v_grid[2], v_grid[2]+rv_range]

                rounds += 1

        ind = np.argsort(v_all)
        v_all = np.asarray(v_all)[ind]
        rms_all = np.asarray(rms_all)[ind]

        # if best RV is too far away from start, just use values around minimum
        # not needed for good frequency resolution of tpl
        pmin = np.argmin(rms_all)
        if abs(v_all[pmin]-rv_guess) > 1:
            rms_all = rms_all[pmin-2:pmin+3]
            v_all = v_all[pmin-2:pmin+3]

        # polyfit through the rms values
        v_gr = np.linspace(v_all[0], v_all[-1], 100)
        pol, resi, _, _, _ = np.polyfit(v_all, rms_all, 2, w=1./np.sqrt(rms_all), full=True)
        sp = np.poly1d(pol)(v_gr)
        rvs = v_gr[np.argmin(sp)]

        gplot.RV2title(", v=%.2f ± %.2f m/s" % (rvs*1000, resi*1000))
        gplot.xlabel('"RV [km/s]"')
        gplot.ylabel('"rms"')
        gplot(v_gr, sp, 'w l lc 9 t "polynomial fit",', v_all, rms_all, 'lc 3 ps 1 pt 2 t "curvefit"')
        # pause()

        # p, e_params = S_mod.fit(pixel_ok, spec_obs_ok, rvs, a, par_wave_guess, s, c=cc, c0=c0, dx=0.1)
        par_rv = rvs
    
    show = (order in look) or (order in lookfast)

    if 1:
        # par from prefit, (not pre-clip)
        par.wave = parguess.wave   # why?
        if ipB:
            par.bkg = [(0, 0)]
            par.ipB = [(ipB[0], 0)]
        if deg_bkg:
            par.bkg = [0]

        par4, e_params = S_mod.fit(pixel_ok, spec_obs_ok, par, dx=0.1*show, sig=sig[i_ok], res=(not createtpl)*show, rel_fac=createtpl*show)
        par = par4

    if kapsig[-1]:
        # second kappa sigma clipping of outliers
        smod = S_mod(pixel, **par)
        resid = spec_obs - smod
        resid[flag_obs != 0] = np.nan

        nr_k1 = np.count_nonzero(flag_obs)
        flag_obs[abs(resid) >= (kapsig[-1]*np.nanstd(resid))] |= flag.clip
        nr_k2 = np.count_nonzero(flag_obs)

        # test if outliers were flagged
        if nr_k1 != nr_k2:
            i_ok = np.where(flag_obs == 0)[0]
            pixel_ok = pixel[i_ok]
            wave_obs_ok = wave_obs[i_ok]
            spec_obs_ok = spec_obs[i_ok]

            par5, e_params = S_mod.fit(pixel_ok, spec_obs_ok, par3, dx=0.1*show, sig=sig[i_ok], res=(not createtpl)*show, rel_fac=createtpl*show)
            par = par5

    if createtpl:
        if tplname:
            # model just the tellurics; exclude stellar lines
            S_star = lambda x: 0*x + 1
            S_mod = model(S_star, lnwave_j, spec_cell_j, specs_molec, IP, **modset)
            
        # modeled telluric spectrum
        spec_model = np.nan * np.empty_like(pixel)
        spec_model[iset] = S_mod(pixel[iset], **par)
        spec_model /= np.nanmedian(spec_model[iset])

        # telluric corrected spectrum
        spec_cor = spec_obs / spec_model
        err_cor = err_obs / spec_model

        # remove regions with strong telluric lines
        spec_cor[spec_model<0.2] = np.nan
        #spec_cor[spec_cor<3*err_cor] = np.nan
        spec_cor[spec_cor<0.01] = np.nan

        wave_model = np.poly1d(par.wave[::-1])(pixel-xcen)
        spec_cor = np.interp(wave_model, wave_model*(1+berv/c)/(1+par.rv/c), spec_cor/np.nanmedian(spec_cor))

        # downweighting by telluric spectrum and errors
        weight = spec_model / (err_cor/np.nanmedian(spec_cor))**2
        # weight[spec_model<0.2] = 0.00001   # downweight deep telluric lines
        weight = np.interp(wave_model, wave_model*(1+berv/c)/(1+par.rv/c), weight)

        # save telluric corrected spectrum
        spec_all[o, 0][n] = wave_model   # updated wavelength
        spec_all[o, 1][n] = spec_cor     # telluric corrected spectrum
        spec_all[o, 2][n] = weight       # weighting for combination of spectra

    if show:
        # overplot flagged and clipped data
        gplot+(pixel[flag_obs != 0], wave_obs[flag_obs != 0], spec_obs[flag_obs != 0], 1*(flag_obs[flag_obs != 0] == flag.clip), 'us (lam?$2:$1):3:(int($4)?5:9) w p pt 6 ps 0.5 lc 9 t "flagged and clipped"')

    if infoprec:
        # estimate velocity precision limit from stellar information content
        # without iodine cell (smoothed to a constant)
        S_pure = model(S_star, lnwave_j, spec_cell_j*0+np.nanmean(spec_cell_j), specs_molec, IP, **modset)
        dS = S_pure(pixel+0.1, **par) - S_pure(pixel, **par)   # flux gradient from finite difference
        du = 1000 * c * np.diff(wave_obs)*0.1 / wave_obs[:-1]   # [m/s] velocity differential from initial solution
        # assuming spectrum given in photon counts (until viper propagates flux uncertainties)
        varS = abs(spec_obs) + 5**2   # (5 = readout noise)
        # RV precision Eq. (6) from Butler+ (1996PASP..108..500B)
        ev_star = np.sum(((dS[:-1]/du)**2 / varS[:-1])[i_ok])**-0.5
        print(f'Stellar RV precision limit: {ev_star} m/s')

        # estimate velocity precision limit for the iodine from its RV information content
        # now the star is smoothed
        tpl_smooth = np.cumsum(spec_tpl[order])
        wz = 1000   # window size
        tpl_smooth = (tpl_smooth[wz:] - tpl_smooth[:-wz]) / wz
        # gplot(spec_tpl[order], ',', tpl_smooth)
        S_smooth = lambda x: np.interp(x, np.log(wave_tpl[order][wz//2:-wz//2])-berv/c, tpl_smooth)
        iod_pure = model(S_smooth, lnwave_j, spec_cell_j, specs_molec, IP, **modset)
        dS = iod_pure(pixel+0.1, **par) - iod_pure(pixel, **par)   # flux gradient from finite difference
        ev_iod = np.sum(((dS[:-1]/du)**2 / varS[:-1])[i_ok])**-0.5
        print(f'Iodine RV precision limit: {ev_iod} m/s')

        # total precision is the squared sum of both
        # in practice it will be worse, since more parameters are modelled
        ev_total = np.sqrt(ev_star**2 + ev_iod**2)
        print(f'Total RV precision limit: {ev_total} m/s')
        if 1:
            gplot2.xlabel('"pixel"')
            gplot2.ylabel('"flux"')
            gplot2(pixel, spec_obs, flag_obs, f' us 1:2:($3>0?9:1) lc var ps 0.5 t "data ({ev_total:.2f} m/s)",', pixel, iod_pure(pixel, **par)+np.nanmean(spec_obs)/2, flag_obs, f' us 1:2:($3>0?9:2) w l lc var t "offset + IP x iod ({ev_iod:.2f} m/s)",', pixel, S_pure(pixel, **par), flag_obs, f' us 1:2:($3>0?9:3) w l lc var t "IP x star ({ev_star:.2f} m/s)"')
            pause()

    # overplot FTS iodine spectrum
    #gplot+(np.exp(lnwave_j), spec_cell_j/spec_cell_j.max()*spec_obs_ok.max(), 'w l lc 9')
    # overplot stellar spectrum
    #gplot+(np.exp(lnwave_j), S_star(lnwave_j)/S_star(lnwave_j).max()*spec_obs_ok.max(), 'w l lc 9')

    rvo, e_rvo = 1000*par.rv, 1000*par.rv.unc   # convert to m/s
    #prms = S_mod.show([params[0], params[1:1+1+deg_norm], params[2+deg_norm:2+deg_norm+1+deg_wave], params[3+deg_norm+deg_wave:]], pixel_ok, spec_obs_ok, dx=0.1)
    # gplot+(wave_tpl[s_s]*(1-berv/c), spec_tpl[s_s]*parguess_norm, 'w lp lc 4 ps 0.5')
    #gplot+(pixel_ok, S_star(np.log(np.poly1d(b[::-1])(pixel_ok))+(v)/c), 'w lp ps 0.5')
    # gplot+(np.exp(S_star.x), S_star.y, 'w lp ps 0.5 lc 7')

    fmod = S_mod(pixel_ok, **par)
    res = spec_obs_ok - fmod
    prms = np.nanstd(res) / np.nanmean(fmod) * 100
    np.savetxt('res.dat', list(zip(pixel_ok, res)), fmt="%s")

    if order in look:
        pause('look %s:'% o, rvo, '+/- %.2f' % e_rvo)  # globals().update(locals())

    if order in lookres:
        gplot2.palette_defined('(0 "blue", 1 "green", 2 "red")')
        gplot2.var(j=1, lab_ddS='"Finite second derivative 2S(i) - S(i+1) - S(i-1)"')
        # shortcut "j" allows to toggle between flux and second derivative
        gplot2.bind('j "j=(j+1) % 2; set xlabel (j==0? \\"S(i)\\" : lab_ddS) ;repl"')
        gplot2.xlabel('lab_ddS')
        gplot2.ylabel('"residuals S_i - S(i)"')
        gplot2.cblabel('"pixel x_i"')
        #gplot(pixel_ok, spec_obs_ok, S_mod(pixel_ok, **par), 2*spec_obs_ok-f[pixel_ok-1]-f[pixel_ok+1], 'us 3+j:($2-$3):1 w p pt 7 palette t ""')
        gplot2(pixel_ok, spec_obs_ok, S_mod(pixel_ok, **par), 2*S_mod(pixel_ok, **par)-S_mod(pixel_ok-1, **par)-S_mod(pixel_ok+1, **par), 'us 3+j:($2-$3):1 w p pt 7 palette t ""')
        pause(f'lookres {o}')

    if order in lookpar:   
        sa = tplname is not None
        sb = sa + deg_norm+1
        ss = sb + deg_wave+1
        # error estimation
        # uncertainty in continuum
        xl = np.log(np.poly1d(par.wave[::-1])(pixel-xcen))
        Cg = np.poly1d(parguess.norm[::-1])(pixel-xcen)      # continuum guess
        Cp = np.poly1d(par.norm[::-1])(pixel-xcen)    # best continuum 
        X = np.vander(xl, deg_norm+1)[:,::-1].T
        e_Cp = np.einsum('ji,jk,ki->i', X, e_params[sa:sb,sa:sb], X)**0.5
        # uncertainty in wavelength solution
        X = np.vander(xl, deg_wave+1)[:,::-1].T
        lam_g = np.poly1d(parguess.wave[::-1])(pixel-xcen)
        lam = np.poly1d(par.wave[::-1])(pixel-xcen)
        e_lam = np.einsum('ji,jk,ki->i', X, e_params[sb:ss,sb:ss], X)**0.5
        e_wavesol = np.sum((e_lam/lam*3e8)**-2)**-0.5

        # compare the wavelength solutions
        #show_model(i, np.poly1d(b[::-1])(i), np.poly1d(par_wave_guess[::-1])(i), res=True)
        gplot.reset()
        gplot.multiplot("layout 2,2")
        gplot.xlabel('"pixel"').ylabel('"k(x2)"')
        gplot.mxtics().mytics()
        gplot(f'[{ibeg}:{iend}][0:]', pixel, Cg, Cp, e_Cp, 'w l lc 9 t "guess",  "" us 1:3 w l lc 3, "" us 1:($3-$4):($3+$4) w filledcurves fill fs transparent solid 0.2 lc 3 t "1{/Symbol s}" ')
        gplot.xlabel('"pixel"').ylabel('"deviation c * ({/Symbol l} / {/Symbol l}_{guess} - 1) [km/s]"')
        gplot(f'[{ibeg}:{iend}]', pixel, (lam/lam_g-1)*c, ((lam-e_lam)/lam_g-1)*c, ((lam+e_lam)/lam_g-1)*c, 'w l lc 3, "" us 1:3:4 w filledcurves fill fs transparent solid 0.2 lc 3 t "1{/Symbol s}"')
        gplot.xlabel('"[km/s]"').ylabel('"contribution"')
        e_s = e_params[ss,ss]**0.5
        gplot(S_mod.vk, S_mod.IP(S_mod.vk, *parguess.ip), ' lc 9 ps 0.5 t "IP_{guess}", ',
              S_mod.vk, S_mod.IP(S_mod.vk, *par.ip),
                        S_mod.IP(S_mod.vk, *[par.ip[0]-e_s, *par.ip[1:]]),
                        S_mod.IP(S_mod.vk, *[par.ip[0]+e_s, *par.ip[1:]]),
              'lc 3 ps 0.5 t "IP", "" us 1:3:4 w filledcurves fill fs transparent solid 0.2 lc 3 t "1{/Symbol s}"')
        gplot.unset('multiplot')
        pause('lookpar', par.ip)

    return rvo, e_rvo, bjd.jd, berv, par, e_params, prms


obsnames = np.array(sorted(glob.glob(obspath)))[nset]
obsnames = [x for x in obsnames if not any(pat in os.path.basename(x) for pat in nexcl)]

N = len(obsnames)
if not N: pause('no files: ', obspath)

if targname:
    targ = Targ(targname, csv=tag+'.targ.csv').sc

orders = np.r_[oset]
print(orders)

rv = np.nan * np.empty(chunks*len(orders))
e_rv = np.nan * np.empty(chunks*len(orders))

rvounit = open(tag+'.rvo.dat', 'w')
parunit = open(tag+'.par.dat', 'w')

# file headers
colnums = orders if chunks == 1 else [f'{order}-{ch}' for order in orders for ch in range(chunks)]

print('BJD RV e_RV BERV', *map("rv{0} e_rv{0}".format, colnums), 'filename', file=rvounit)

# estimate wavelength range from observation
pixel, wave0, spec0, err0, flag0, bjd, berv = Spectrum(obsnames[0], order=orders[0])
pixel, wave1, spec1, err1, flag1, bjd, berv = Spectrum(obsnames[0], order=orders[-1])
    
obs_lmin = np.min([wave0[0], wave0[-1], wave1[0], wave1[-1]])
obs_lmax = np.max([wave0[0], wave0[-1], wave1[0], wave1[-1]])

####  FTS  ####
# using the supersampled log(wavelength) space with knot index j

if ftsname != 'None':
    wave_cell, spec_cell, lnwave_j_full, spec_cell_j_full = FTS(ftsname)
else:
    # create fake cell spectrum 
    wave_cell = np.linspace(obs_lmin, obs_lmax, len(pixel)*len(orders)*100)
    spec_cell = wave_cell*0 + 1
    u = np.log(wave_cell)
    lnwave_j_full = np.arange(u[0], u[-1], 100/3e8)
    spec_cell_j_full = lnwave_j_full*0 + 1 

if nocell:
    # option nocell will be removed in near future
    spec_cell = spec_cell*0 + 1
    spec_cell_j_full = spec_cell_j_full*0 + 1

# mask wavelengths with strong tellurics
mskatm = lambda x: np.interp(x, *np.genfromtxt(viperdir+'lib/mask_vis1.0.dat').T)

if flagfile:
    # user created file for removal of selected regions
    msk = np.genfromtxt(flagfile, names=True, invalid_raise=False, missing_values = {'order':"-"}, filling_values={'order':np.nan}, delimiter=' ').view(np.recarray)
    
    msk_o = msk[np.isfinite(msk.order)]		# selected pixel in order
    msk_l = msk[np.isnan(msk.order)]		# selected wavelength/lambda ranges

    if len(msk_l):
        msk_w = np.asarray(np.concatenate([msk_l.start, msk_l.end, msk_l.start-0.05, msk_l.end+0.05]))
        msk_f = np.asarray(np.concatenate([np.ones(2*len(msk_l.start)), np.zeros(2*len(msk_l.start))]))
    
        ind = np.argsort(msk_w)
        msk_w, msk_f = msk_w[ind], msk_f[ind]
        

#### Telluric model ####
if 'add' in telluric:
    # read in telluric spectra for wavelength range of the instrument

    bands_all = ['vis', 'J', 'H', 'K']
    wave_band = [0, 9000, 14000, 18500]
    
    # select which bands are covered by the observation
    w0 = obs_lmin - wave_band
    w1 = obs_lmax - wave_band
    bands = bands_all[np.argmin(w0[w0 >= 0]): int(np.argmin(w1[w1 >= 0]) + 1)]

    specs_molec_all = defaultdict(list)
    wave_atm_all = defaultdict(list)
    molec_sel = molec
    
    for band in bands:

        hdu = fits.open(viperdir+'/lib/atmos/stdAtmos_'+band+'.fits')
        cols = hdu[1].columns.names
        data = hdu[1].data

        if molec_sel[0] == 'all': molec = cols[1:]
        
        # add wavelength shift
        # synthetic telluric spectra (molecfit) are laboratory wavelengths
        # shift was determined empirical from several observations                    
        for i_mol, mol in enumerate(molec):
            if (mol != 'lambda') and (mol in cols): 
                specs_molec_all[mol].extend(data[mol])
                wave_atm_all[mol].extend(data['lambda'] * (1 + (-0.249/3e5)))
                
        molec = np.array(list(specs_molec_all.keys()))

# collect all spectra for createtpl function
spec_all = defaultdict(dict)

####  stellar template  ####

if tplname:
    print('reading stellar template')
    wave_tpl, spec_tpl = {}, {}
    for order in orders:
        wave_tplo, spec_tplo = Tpl(tplname, order=order, targ=targ)
        if oversampling:
            us = np.linspace(np.log(wave_tplo[0]), np.log(wave_tplo[-1]), oversampling*wave_tplo.size)
            spec_tplo = np.nan_to_num(spec_tplo)
            fs = CubicSpline(np.log(wave_tplo), spec_tplo)(us)
            wave_tpl[order], spec_tpl[order] = np.exp(us), fs
        else:
            wave_tpl[order], spec_tpl[order] = wave_tplo, spec_tplo
else:
    # no template given; model pure iodine
    wave_tpl, spec_tpl = [wave_cell[[0, -1]]]*100, [np.ones(2)]*100


if createtpl:
    wave_tplo, spec_tplo = Tpl(obsnames[-1], order=orders[-1], targ=targ)
    wmax = np.max(wave_tplo)
else:
    wmax = np.max(wave_tpl[orders[-1]])

if telluric == 'add' and (wave_cell[-1] < wmax):
    # extend wavelength range for telluric modelling
    # iodine ends around order 36 for TLS and OES
    # at higher orders modelling with telluric lines instead of iodine is possible

    wave_cell_ext = np.arange(wave_cell[-1], wmax, wave_cell[-1]-wave_cell[-2])[1:]
    spec_cell_ext = np.ones_like(wave_cell_ext)

    wave_cell = np.append(wave_cell, wave_cell_ext)
    spec_cell = np.append(spec_cell, spec_cell_ext)

    lnwave_j = np.log(wave_cell)
    lnwave_j_full = np.arange(lnwave_j[0], lnwave_j[-1], 100/3e8)
    spec_cell_j_full = np.interp(lnwave_j_full, lnwave_j, spec_cell)

    if not tplname:
        wave_tpl, spec_tpl = [wave_cell[[0, -1]]]*100, [np.ones(2)]*100

T = time.time()
headrow = True
for n, obsname in enumerate(obsnames):
    if n == 0:
        # clear up the residual directory
        if os.path.isdir(viperdir+'res') and os.listdir(viperdir+'res'):
            os.system('rm -rf '+viperdir+'res/*.dat')
    filename = os.path.basename(obsname)
    print(f"{n+1:3d}/{N}", filename)
    for i_o, o in enumerate(orders):
        for ch in np.arange(chunks):
            try:
                gplot.RV2title = lambda x: gplot.key('title noenhanced "%s (n=%s, o=%s%s)"'% (filename, n+1, o, x))
                gplot.RV2title('')
        
                rv[i_o*chunks+ch], e_rv[i_o*chunks+ch], bjd, berv, params, e_params, prms = fit_chunk(o, ch, obsname=obsname, targ=targ)

                print(n+1, o, ch, rv[i_o*chunks+ch], e_rv[i_o*chunks+ch])
                # just for compability, remove Params(ipB=[]) later !!
                if 'ipB' in params: params.pop('ipB')
                if not deg_bkg: params.pop('bkg', None)                       
                params.rv.value *= 1000.   # convert to m/s -> same unit in .par.dat and .rvo.dat     
                params.rv.unc *= 1000.

                if headrow:
                    headrow = False
                    colnames = ["".join(map(str,x)) for x in params.flat().keys()]
                    print('BJD n order chunk', *map("{0} e_{0}".format, colnames), 'prms', file=parunit)

                flat_params = [f"{d.value} {d.unc}" for d in params.flat().values()]
                print(bjd, n+1, o, ch, *flat_params, prms, file=parunit)
                # store residuals
                os.system('mkdir -p res; touch res.dat')
                os.system('mv res.dat res/%03d_%03d.dat' % (n, o))

            except Exception as e:
                if repr(e) == 'BdbQuit()':
                    exit()
                print("Order failed due to:", repr(e))

    if not np.isnan(rv).all():
        oo = np.isfinite(e_rv)
        if oo.sum() == 1:
            RV = rv[oo][0]
            e_RV = e_rv[oo][0]
        else:
            RV = np.nanmean(rv[oo])
            e_RV = np.nanstd(rv[oo])/(oo.sum()-1)**0.5
        print('RV:', RV, e_RV, bjd, berv)

        print(bjd, RV, e_RV, berv, *sum(zip(rv, e_rv), ()), filename, file=rvounit)
        print(file=parunit)


if createtpl:
    # combine all telluric corrected spectra to a final template
    wave_tpl_new = {}
    spec_tpl_new = {}
    err_tpl_new = {}
    orders_ok = sorted(set([kk[0] for kk in spec_all.keys()]))
    for order in orders_ok:
        gplot.reset()
        gplot.key("title 'order: %s' noenhance" % (order))
        gplot.xlabel('"Vacuum wavelength [Å]"')
        gplot.ylabel('"flux"')
        gplot.yrange("[%g:%g]" % (-1, 1.6))
        wave_t = np.array(list(spec_all[order,0].values()))     # wavelength
        spec_t = np.array(list(spec_all[order,1].values()))     # data
        weight_t = np.array(list(spec_all[order,2].values()))   # weighting
        weight_t[np.isnan(spec_t)] = 0
        weight_t[spec_t<0] = 0
        # weight_t[spec_t>1.15] = np.nanmin(weight_t)/10.

      #  spec_tpl_new[order] = np.nansum(spec_t*weight_t, axis=0) / np.nansum(weight_t, axis=0)
        #wave_tpl_new[order] = np.nanmean(wave_t, axis=0)
        wave_tpl_new[order] = wave_t[0]

        if len(spec_t) > 1:
            # combine several observations to one tpl
            for nn in range(1, len(spec_t)):
                # flag outlier points and spectra
                valid = np.isfinite(spec_t[nn])
                spec_cubic = CubicSpline(wave_t[nn][valid], spec_t[nn][valid])(wave_t[0])
                spec_cubic[valid==0] = np.nan
                spec_t[nn] = spec_cubic

                # weight_cubic = CubicSpline(wave_t[nn][valid], weight_t[nn][valid])(wave_t[0])
                # weight_t[nn][valid] = weight_cubic[valid]
                weight_t[nn] = np.interp(wave_t[0], wave_t[nn][valid], weight_t[nn][valid])
                weight_t[nn][valid==0] = np.nan
                weight_t[nn][weight_t[nn]==0] = np.nan

            if kapsig_ctpl:
                spec_mean = np.nansum(spec_t*weight_t, axis=0) / np.nansum(weight_t, axis=0)
                for nn in range(0, len(spec_t)):
                    weight_t[nn][np.abs(spec_t[nn]-spec_mean)>kapsig_ctpl] = np.nan

            spec_tpl_new[order] = np.nansum(spec_t*weight_t, axis=0) / np.nansum(weight_t, axis=0)
            err_tpl_new[order] = np.nanstd(spec_t, axis=0)
          
        else:
            spec_tpl_new[order] = spec_t[0]
            err_tpl_new[order] = spec_t[0]*np.nan

        if (order in lookfast) or (order in look) or (order in lookctpl):
            gplot(wave_tpl_new[order], spec_tpl_new[order] - 1 , 'w l lc 7 t "combined tpl"')
            for n in range(len(spec_t)):
                gplot+(wave_tpl_new[order], spec_t[n]/np.nanmedian(spec_t[n]), 'w l t "%s"' % (os.path.split(obsnames[n])[1]))          
            #gplot+(wave_tpl_new[order], np.nanstd(spec_t, axis=0)+1.5, 'w l t ""')
        if (order in look) or (order in lookctpl):
            pause()

    Inst.write_fits(wave_tpl_new, spec_tpl_new, err_tpl_new, obsnames, tag)

rvounit.close()
parunit.close()
convert_output.convert_data(tag, args, dat='dat' in oformat, fits='fits' in oformat, cpl='cpl' in oformat)

T = time.time() - T
Tfmt = lambda t: time.strftime("%Hh%Mm%Ss", time.gmtime(t))
print("processing time total:       ", Tfmt(T))
print("processing time per spectrum:", Tfmt(T/N))
print("processing time per chunk:   ", Tfmt(T/N/orders.size))

if 'cpl' in oformat or 'fits' in oformat:
    tag += '_rvo_par.fits'
else:
    tag += '.rvo.dat'

if not createtpl:
    vpr = vpr.VPR(tag)   # to print info statistic
    if len(lookfast) or len(look):
        gplot.reset()
        vpr.plot_RV()
print(tag, 'done.')

def run():
    pass
