import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf
from astropy.modeling.models import Voigt1D

from gplot import *

c = 3e5   # [km/s] speed of light

# IP sampling in velocity space
# index k for IP space
def IP(vk, s=2.2):
    """Gaussian IP"""
    IP_k = np.exp(-(vk/s)**2/2)
    IP_k /= IP_k.sum()          # normalise IP
    return IP_k

def IP_sg(vk, s=2.2, e=2.):
    """super Gaussian"""
    IP_k = np.exp(-abs(vk/s)**e)
    IP_k /= IP_k.sum()          # normalise IP
    return IP_k

def IP_ag(vk, s=2.2, a=0):
    '''
    Asymmetric (skewed) Gaussian.

    Example
    -------
    >>> vk = np.arange(-50, 50+1)
    >>> gplot(vk, IP_ag(vk, s=10.), IP_ag(vk, s=10., a=100), ', "" us 1:3, "" us 1:($3-$2)*0.5, 0')

    '''
    b = a / np.sqrt(1+a**2) * np.sqrt(2/np.pi)
    ss = s / np.sqrt(1-b**2)    # readjust scale parameter to have same variance as unskewed Gaussian
    vk = (vk + ss*b) / ss       # recenter to have zero mean
    IP_k = np.exp(-vk**2/2) * (1+erf(a/np.sqrt(2)*vk))   # Gauss IP * erf
    IP_k /= IP_k.sum()          # normalise IP
    return IP_k

def IP_mg(vk, s0=2, a1=0.1):
    """IP for multiple, zero-centered Gaussians."""
    s1 = 4 * s0   # width of second Gaussian with fixed relative width
    a1 = a1 / 10  # relative ampitude
    IP_k = np.exp(-(vk/s0)**2)         # main Gaussian
    IP_k += a1 * np.exp(-(vk/s1)**2)   # a satellite Gaussian
    IP_k = IP_k.clip(0, None)
    IP_k /= IP_k.sum()          # normalise IP
    return IP_k


IPs = {'g': IP, 'sg': IP_sg, 'ag': IP_ag, 'mg': IP_mg, 'bnd': 'bnd'}


def poly(x, a):
    # redefine polynomial for argument order and adjacent coefficients
    return np.polyval(a[::-1], x)

def pade(x, a, b):
    '''
    rational polynomial
    b: denominator coefficients b1, b2, ... (b0 is fixed to 1)
    '''
    y = poly(x, a) / (1+x*poly(x, b))
    return y


class model:
    '''
    The forward model.

    '''
    def __init__(self, *args, envelope=poly, IP_hs=50, icen=0):
        # IP_hs: Half size of the IP (number of sampling knots).
        # icen : Central pixel (to center polynomial for numeric reason).

	self.icen = icen
        self.S_star, self.lnwave_j, self.spec_cell_j, self.spec_atm, self.IP = args
        # convolving with IP will reduce the valid wavelength range
        self.dx = self.lnwave_j[1] - self.lnwave_j[0]   # step size of the uniform sampled grid
        self.IP_hs = IP_hs
        self.vk = np.arange(-IP_hs, IP_hs+1) * self.dx * c
        self.lnwave_j_eff = self.lnwave_j[IP_hs:-IP_hs]    # valid grid
        self.envelope = envelope
        #print("sampling [km/s]:", self.dx*c)

    def __call__(self, pixel, rv, coeff_norm, coeff_wave, coeff_ip, coeff_atm, coeff_bkg=[0]):
        # wavelength solution
        #    lam(x) = b0 + b1 * x + b2 * x^2
        lnwave_obs = np.log(np.poly1d(coeff_wave[::-1])(pixel-self.icen))

        # IP convolution
        if len(self.spec_atm) == 0:
            Sj_eff = np.convolve(self.IP(self.vk, *coeff_ip), self.S_star(self.lnwave_j-rv/c) * (self.spec_cell_j + coeff_bkg[0]), mode='valid')
        else:
            # telluric forward modelling
            atm = np.ones(len(self.spec_atm[0]))
            for coeff_mol, sp in zip(coeff_atm, self.spec_atm):
                atm *= sp**np.abs(coeff_mol)

	    # variable telluric wavelength shift; one shift for all molecules
            if len(coeff_atm) == len(self.spec_atm)+1:
                atm = np.interp(self.lnwave_j, self.lnwave_j-np.log(1+coeff_atm[-1]/c), atm)

            Sj_eff = np.convolve(self.IP(self.vk, *coeff_ip), self.S_star(self.lnwave_j-rv/c) * (self.spec_cell_j * atm + coeff_bkg[0]), mode='valid')

        # sampling to pixel
        Si_eff = np.interp(lnwave_obs, self.lnwave_j_eff, Sj_eff)

        # flux normalisation
        Si_mod = self.envelope(pixel-self.icen, coeff_norm) * Si_eff
        #Si_mod = self.envelope((np.exp(lnwave_obs)-b[0]-coeff_norm[-1]), coeff_norm[:-1]) * Si_eff
        return Si_mod

    def fit(self, pixel, spec_obs, par_rv=None, par_norm=[], par_wave=[], par_ip=[], par_atm=[], par_bkg=[], parfix_rv=None, parfix_norm=[], parfix_wave=[], parfix_ip=[], parfix_atm=[], parfix_bkg=[], sig=[], **kwargs):
        '''
        Generic fit wrapper.
        '''
        parfix_rv = () if parfix_rv is None else (parfix_rv,)
        par_rv = () if par_rv is None else (par_rv,)
        sv = slice(0, len(par_rv))
        sa = slice(sv.stop, sv.stop+len(par_norm))
        sb = slice(sa.stop, sa.stop+len(par_wave))
        ss = slice(sb.stop, sb.stop+len(par_ip))
        st = slice(ss.stop, ss.stop+len(par_atm))
        sc = slice(st.stop, st.stop+len(par_bkg))
        parfix_norm = tuple(parfix_norm)
        parfix_wave = tuple(parfix_wave)
        parfix_ip = tuple(parfix_ip)
        parfix_atm = tuple(parfix_atm)
        parfix_bkg = tuple(parfix_bkg)

        S_model = lambda x, *params: self(x, *params[sv]+parfix_rv, params[sa]+parfix_norm, params[sb]+parfix_wave, params[ss]+parfix_ip, params[st]+parfix_atm, params[sc]+parfix_bkg)

        params, e_params = curve_fit(S_model, pixel, spec_obs, p0=[*par_rv, *par_norm, *par_wave, *par_ip, *par_atm, *par_bkg], sigma=sig, absolute_sigma=False, epsfcn=1e-12)

        par_rv = (*params[sv], *np.diag(e_params)[sv])
        params = tuple(params)
        params = [*params[sv]+parfix_rv, params[sa]+parfix_norm, params[sb]+parfix_wave, params[ss]+parfix_ip, params[st]+parfix_atm, params[sc]+parfix_bkg]

        if kwargs:
            self.show(params, pixel, spec_obs, par_rv=par_rv, **kwargs)
        return params, e_params

    def show(self, params, x, y, par_rv=None, res=True, x2=None, dx=None, rel_fac=None):
        '''
        res: Show residuals.
        x2: Values for second x axis.
        rel_fac: Factor for relative residuals.
        dx: Subpixel step size for the model [pixel].
        '''
        ymod = self(x, *params)
        if x2 is None:
            x2 = np.poly1d(params[2][::-1])(x-self.icen)
        if par_rv:
            gplot.RV2title(", v=%.2f ± %.2f m/s" % (par_rv[0]*1000, np.sqrt(par_rv[1])*1000))

        gplot.put("if (!exists('lam')) {lam=1}")

        gplot.key('horizontal')
        gplot.xlabel('lam?"Vacuum wavelength [Å]":"Pixel x"')
        gplot.ylabel('"flux"')
        # toggle between pixel and wavelength with shortcut "$"
        gplot.bind('"$" "lam=!lam; set xlabel lam?\\"Vacuum wavelength [Å]\\":\\"Pixel x\\"; replot"')
        args = (x, y, ymod, x2, 'us lam?4:1:2:3 w lp pt 7 ps 0.5 t "S_i",',
          '"" us lam?4:1:3 w p pt 6 ps 0.5 lc 3 t "S(i)"')
        prms = np.nan   # percentage prms
        if dx:
            xx = np.arange(x.min(), x.max(), dx)
            xx2 = np.poly1d(params[2][::-1])(xx-self.icen)
            yymod = self(xx, *params)
            args += (",", xx, yymod, xx2, 'us lam?3:1:2 w l lc 3 t ""')
        if res or rel_fac:
            # linear or relative residuals
            col2 = rel_fac * np.mean(ymod) * (y/ymod - 1) if rel_fac else y - ymod
            rms = np.std(col2)
            prms = rms / np.mean(ymod) * 100
            gplot.mxtics().mytics().my2tics()
        if res:
            args += (",", x, col2, x2, "us lam?3:1:2 w p pt 7 ps 0.5 lc 1 t 'res (%.3g \~ %.3g%%)', 0 lc 3 t ''" % (rms, prms))
        if rel_fac:
            args += (",", x, col2, x2, "us lam?3:1:2 w l lc 1 t 'res (%.3g \~ %.3g%%)', 0 lc 3 t ''" % (rms, prms))

        gplot(*args)
        return prms


class model_bnd(model):
    '''
    The forward model with band matrix.

    '''
    def base(self, x=None, degk=3, sig_k=1):
        '''
        Setup the base functions.
        '''
        self.x = x
        bg = self.IP
        # the initial trace
        lnwave_obs = np.log(np.poly1d(bg[::-1])(x-self.icen))
        j = np.arange(self.lnwave_j.size)
        jx = np.interp(lnwave_obs, self.lnwave_j, j)
        # position of Gaussians
        vl = np.array([-1, 0, 1])[np.newaxis,np.newaxis,:]
        vl = np.array([-1.4, -0.7, 0, 0.7, 1.4])[np.newaxis,np.newaxis,:]
        #vl = np.array([0])[np.newaxis,np.newaxis,:]

        # bnd -- lists all j's contributing to x
        self.bnd = jx[:,np.newaxis].astype(int) + np.arange(-self.IP_hs, self.IP_hs+1)

        # base for multi-Gaussians
        self.BBxjl = np.exp(-(self.lnwave_j[self.bnd][...,np.newaxis]-lnwave_obs[:,np.newaxis,np.newaxis]+sig_k*vl)**2/sig_k**2)

        # base function for flux polynomial
        self.Bxk = np.vander(x-x.mean(), degk)[:,::-1]
        #Bxk = np.vander(jx-jx.mean(), len(ak))[:,::-1]   # does not provide proper output for len(ak) > 4

    def Axk(self, v, **kwargs):
        if kwargs:
            self.base(**kwargs)
        starj = self.S_star(self.lnwave_j-v/c)  # np.interp(self.lnwave_j, np.log(tpl_w)-berv/3e5, tpl_f/np.median(tpl_f))
        _Axkl = np.einsum('xj,xjl,xk->xkl', (starj*self.spec_cell_j)[self.bnd], self.BBxjl, self.Bxk)
        return _Axkl

    def IPxj(self, akl, **kwargs):
        if kwargs:
            self.base(**kwargs)
        # sum all Gaussians
        IPxj = np.einsum('xjl,xk,kl->xj', self.BBxjl, self.Bxk, akl.reshape(self.Bxk.shape[1], -1)) #.reshape(AAxkl[0].shape))
        return IPxj

    def fit(self, f, v, **kwargs):
        Axkl = self.Axk(v, **kwargs)
        return np.linalg.lstsq(Axkl.reshape((len(Axkl), -1)), f, rcond=1e-32)

    def __call__(self, x, v, ak, **kwargs):
        Axkl = self.Axk(v, **kwargs)
        # dl = np.array([0.5, 2, 0.5])   # amplitude of Gaussian
        #ak = np.array([1, 0, 0, 0, 0])
        #fx = Axk @ ak
        # fx = AAxkl.reshape((len(AAxkl), -1)) @ aakl
        fx = Axkl.reshape((len(Axkl), -1)) @ ak
        return fx


def show_model(x, y, ymod, res=True):
    gplot(x, y, ymod, 'w lp pt 7 ps 0.5 t "S_i",',
          '"" us 1:3 w lp pt 6 ps 0.5 lc 3 t "S(i)"')
    if res:
        rms = np.std(y-ymod)
        gplot.mxtics().mytics().my2tics()
        # overplot residuals
        gplot.y2range('[-0.2:2]').ytics('nomirr').y2tics()
        gplot+(x, y-ymod, "w p pt 7 ps 0.5 lc 1 axis x1y2 t 'res %.3g', 0 lc 3 axis x1y2" % rms)
