import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf
from astropy.modeling.models import Voigt1D

from gplot import *

c = 3e5   # [km/s] speed of light

# IP sampling in velocity space
# index k for IP space
def IP(vk, s=2.2):
    IP_k = np.exp(-(vk/s)**2/2)   # Gauss IP
    #IP_k += 0.07*np.exp(-((vk+1.)/s)**2)   # an asymmetry Gauss IP
    IP_k /= IP_k.sum()          # normalise IP
    return IP_k

def IP_sg(vk, s=2.2, e=2.):
    # super Gaussian
    IP_k = np.exp(-abs(vk/s)**e)   # Gauss IP
    IP_k /= IP_k.sum()          # normalise IP
    return IP_k

def IP_ag(vk, s=2.2, a=0):
    '''
    Asymmetric (skewed) Gaussian.
    
    Example
    -------
    >>> vk = np.arange(-50, 50+1) 
    >>> gplot(vk, IP_ag(vk, s=10.), IP_ag(vk, s=10., a=100),', "" us 1:3, "" us 1:($3-$2)*0.5, 0')

    '''
    b = a / np.sqrt(1+a**2) * np.sqrt(2/np.pi)
    ss = s / np.sqrt(1-b**2)    # readjust scale parameter to have same variance as unskewed Gaussian
    vk = (vk + ss*b) / ss       # recenter to have zero mean
    IP_k = np.exp(-vk**2/2) * (1+erf(a/np.sqrt(2)*vk))  # Gauss IP * erf
    IP_k /= IP_k.sum()          # normalise IP
    return IP_k


def IP_mg(vk, s0=2, a1=0.1):
    ''' IP for multiple, zero-centered Gaussians. '''
    #print(s)
    s1 = 4 * s0   # width of second Gaussian with fixed relative width
    a1 = a1 /10   # relative ampitude
    IP_k = np.exp(-(vk/s0)**2)   # Gauss IP
    IP_k += a1*np.exp(-(vk/s1)**2)   # Gauss IP
    IP_k = IP_k.clip(0,None)
    IP_k /= IP_k.sum()          # normalise IP
    return IP_k

IPs = {'g':IP, 'sg': IP_sg, 'ag': IP_ag, 'mg':IP_mg, 'bnd': 'bnd'}

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
    The forward model

    '''
    def __init__(self, *args, envelope=poly, IP_hs=50, icen=0):
        # IP_hs: Half size of the IP (number of sampling knots).
        # icen : Central pixel (to center polynomial for numeric reason).
        self.icen = icen
        self.S_star, self.uj, self.iod_j, self.atmj, self.IP = args
        # convolving with IP will reduce the valid wavelength range
        self.dx = self.uj[1] - self.uj[0]  # sampling in uniform resampled Iod
        self.IP_hs = IP_hs
        self.vk = np.arange(-IP_hs, IP_hs+1) * self.dx * c
        self.uj_eff = self.uj[IP_hs:-IP_hs]
        self.envelope = envelope
        #print("sampling [km/s]:", self.dx*c)

    def __call__(self, i, v, a, b, s,t, cc=[0]):
        # wavelength solution 
        #    lam(x) = b0 + b1 * x + b2 * x^2
        ui = np.log(np.poly1d(b[::-1])(i-self.icen))

        # IP convolution
        if len(self.atmj) == 0:
            Sj_eff = np.convolve(self.IP(self.vk, *s), self.S_star(self.uj-v/c) * (self.iod_j + cc[0]), mode='valid')
        else:
            # telluric forward modelling
            atm = np.ones(len(self.atmj[0]))
            for aj in range(0,len(self.atmj),1):
                atm *= (self.atmj[aj]**np.abs(t[aj]))
            Sj_eff = np.convolve(self.IP(self.vk, *s), self.S_star(self.uj-v/c) * (self.iod_j * atm + cc[0]), mode='valid')

        # sampling to pixel
        Si_eff = np.interp(ui, self.uj_eff, Sj_eff)

        # flux normalisation
        Si_mod = self.envelope(i-self.icen, a) * Si_eff
        #Si_mod = self.envelope((np.exp(ui)-b[0]-a[-1]), a[:-1]) * Si_eff
        return Si_mod
 
    def fit(self, x, f, v=None, a=[], b=[], s=[],t=[], c=[], v0=None, a0=[], b0=[], s0=[],t0=[], c0=[],sig=[], **kwargs):
        '''
        Generic fit wrapper.
        '''
        v0 = () if v0 is None else (v0,)
        v = () if v is None else (v,)
        sv = slice(0, len(v))
        sa = slice(sv.stop, sv.stop+len(a))
        sb = slice(sa.stop, sa.stop+len(b))
        ss = slice(sb.stop, sb.stop+len(s))
        st = slice(ss.stop, ss.stop+len(t))
        sc = slice(st.stop, st.stop+len(c))
        a0 = tuple(a0)
        b0 = tuple(b0)
        s0 = tuple(s0)
        t0 = tuple(t0)
        c0 = tuple(c0)
        S = lambda x, *p: self(x, *p[ sv]+v0, p[sa]+a0, p[sb]+b0, p[ss]+s0, p[st]+t0, p[sc]+c0)
        p, e_p = curve_fit(S, x, f, p0=[*v, *a, *b, *s,*t, *c],sigma=sig, absolute_sigma=False, epsfcn=1e-12)
        v = (*p[sv], *np.diag(e_p)[sv])
        p = tuple(p)
        p = [*p[sv]+v0, p[sa]+a0, p[sb]+b0, p[ss]+s0, p[st]+t0, p[sc]+c0]
        if kwargs:
            self.show(p, x, f, v=v, **kwargs)
        return p, e_p

    def show(self, p, x, y, v=None, res=True, x2=None, dx=None):
        '''
        res: Show residuals.
        x2: Values for second x axis.
        dx: Subpixel step size for the model [pixel].

        '''
        ymod = self(x, *p)
        if x2 is None:
            x2 = np.poly1d(p[2][::-1])(x-self.icen)
        if v:
            gplot.RV2title(", v=%.2f ± %.2f m/s" % (v[0]*1000, np.sqrt(v[1])*1000))

        gplot.put("if (!exists('lam')) {lam=1}")

        gplot.key('horizontal')
        gplot.xlabel('lam?"Vaccum wavelength [Å]":"Pixel x"')
        gplot.ylabel('"flux"')
        # toggle between pixel and wavelength with shortcut "$"
        gplot.bind('"$" "lam=!lam; set xlabel lam?\\"Vaccum wavelength [Å]\\":\\"Pixel x\\"; replot"')
        args = (x, y, ymod, x2, 'us lam?4:1:2:3 w lp pt 7 ps 0.5 t "S_i",',
          '"" us lam?4:1:3 w p pt 6 ps 0.5 lc 3 t "S(i)"')
        prms = np.nan   # percentage prms
        if dx:
            xx = np.arange(x.min(), x.max(), dx)
            xx2 = np.poly1d(p[2][::-1])(xx-self.icen)
            yymod = self(xx, *p)
            args += (",", xx, yymod, xx2, 'us lam?3:1:2 w l lc 3 t ""')
        if res:
            rms = np.std(y-ymod)
            prms = rms / np.mean(ymod) * 100
            gplot.mxtics().mytics().my2tics()
            # overplot residuals
            #gplot.y2range('[-0.2:2]').ytics('nomirr').y2tics()
            #args += (",", x, y-ymod, x2, "us lam?3:1:2 w p pt 7 ps 0.5 lc 1 axis x1y2 t 'res (%.3g \~ %.3g%%)', 0 lc 3 axis x1y2 t ''" % (rms, rms/np.mean(ymod)*100))
            args += (",", x, y-ymod, x2, "us lam?3:1:2 w p pt 7 ps 0.5 lc 1 t 'res (%.3g \~ %.3g%%)', 0 lc 3 t ''" % (rms, rms/np.mean(ymod)*100))
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
        ui = np.log(np.poly1d(bg[::-1])(x-self.icen))
        j = np.arange(self.uj.size)
        jx = np.interp(ui, self.uj, j)
        # position of Gaussians
        vl = np.array([-1,0,1])[np.newaxis,np.newaxis,:]
        vl = np.array([-1.4,-0.7,0,0.7,1.4])[np.newaxis,np.newaxis,:]
        #vl = np.array([0])[np.newaxis,np.newaxis,:]

        # bnd -- lists all j's contributing to x
        self.bnd = jx[:,np.newaxis].astype(int) + np.arange(-self.IP_hs, self.IP_hs+1)

        # base for multi-Gaussians
        self.BBxjl = np.exp(-(self.uj[self.bnd][...,np.newaxis]-ui[:,np.newaxis,np.newaxis]+sig_k*vl)**2/sig_k**2)

        # base function for flux polynomial
        self.Bxk = np.vander(x-x.mean(), degk)[:,::-1]
        #Bxk = np.vander(jx-jx.mean(), len(ak))[:,::-1]   # does not provide proper output for len(ak) > 4

    def Axk(self, v, **kwargs):
        if kwargs:
            self.base(**kwargs)
        starj = self.S_star(self.uj-v/c)  # np.interp(self.uj, np.log(tpl_w)-berv/3e5, tpl_f/np.median(tpl_f))
        _Axkl = np.einsum('xj,xjl,xk->xkl', (starj*self.iod_j)[self.bnd], self.BBxjl, self.Bxk)
        return _Axkl

    def IPxj(self, akl, **kwargs):
        if kwargs:
            self.base(**kwargs)
        # sum all Gaussians
        IPxj = np.einsum('xjl,xk,kl->xj', self.BBxjl, self.Bxk, akl.reshape(self.Bxk.shape[1],-1)) #.reshape(AAxkl[0].shape))
        return IPxj

    def fit(self, f, v, **kwargs):
        Axkl = self.Axk(v, **kwargs)
        return np.linalg.lstsq(Axkl.reshape((len(Axkl),-1)), f, rcond=1e-32)

    def __call__(self, x, v, ak, **kwargs):
        Axkl = self.Axk(v, **kwargs)
#        dl = np.array([0.5,2,0.5])   # amplitude of Gaussian
        #ak = np.array([1, 0, 0, 0, 0])
        #fx = Axk @ ak
#        fx = AAxkl.reshape((len(AAxkl),-1)) @ aakl
        fx = Axkl.reshape((len(Axkl),-1)) @ ak
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
