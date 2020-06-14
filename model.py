import numpy as np
from scipy import interpolate

from gplot import *
gplot.tmp = '$'

c = 3e5   # [km/s] speed of light

# IP sampling in velocity space
# index k for IP space
def IP(vk, s=2.2):
    IP_k = np.exp(-(vk/s)**2)   # Gauss IP
    IP_k /= IP_k.sum()          # normalise IP
    return IP_k

def IP2(vk, s):
    ''' IP for multiple, zero-centered Gaussians '''
    print(s)
    s0, a1, s1 = s
    IP_k = np.exp(-(vk/s0)**2)   # Gauss IP
    IP_k += a1*np.exp(-(vk/s1)**2)   # Gauss IP
    IP_k /= IP_k.sum()          # normalise IP
    return IP_k


class model:
    '''
    The forward model

    '''
    def __init__(self, *args, IP_hw=50):
        self.S_star, self.xj, self.iod_j, self.IP = args
        # convolving with IP will reduce the valid wavelength range
        self.dx = self.xj[1] - self.xj[0]  # sampling in uniform resampled Iod
        self.vk = np.arange(-IP_hw,IP_hw+1) * self.dx * c
        self.xj_eff = self.xj[IP_hw:-IP_hw]
        #print("sampling [km/s]:", self.dx*c)

    def __call__(self, i, v, a, b, s):
        # wavelength solution
        xi = np.log(np.poly1d(b[::-1])(i))
        # IP convolution
        Sj_eff = np.convolve(self.IP(self.vk, s), self.S_star(self.xj-v/c) * self.iod_j, mode='valid')
        # sampling to pixel
        Si_eff = interpolate.interp1d(self.xj_eff, Sj_eff)(xi)
        # flux normalisation
        Si_mod = np.poly1d(a[::-1])(xi) * Si_eff
        return Si_mod

    def show(self, p, x, y, res=True, x2=None, dx=None):
        '''
        res: Show residuals.
        x2: Values for second x axis.
        dx: Subpixel step size for the model [pixel].

        '''
        ymod = self(x, *p)
        x2 = np.poly1d(p[2][::-1])(x)
        gplot.put("if (!exists('lam')) {lam=1}")
        
        #gplot.key('horizontal')
        gplot.xlabel('lam?"Vaccum wavelength [Å]":"Pixel x"')
        gplot.ylabel('"flux"')
        # toggle between pixel and wavelength with shortcut "$"
        gplot.bind('"$" "lam=!lam; set xlabel lam?\\"Vaccum wavelength [Å]\\":\\"Pixel x\\"; replot"')
        args = (x, y, ymod, x2, 'us lam?4:1:2:3 w lp pt 7 ps 0.5 t "S_i",',
          '"" us lam?4:1:3 w p pt 6 ps 0.5 lc 3 t "S(i)"')
        if dx:
            xx = np.arange(x.min(), x.max(), dx)
            xx2 = np.poly1d(p[2][::-1])(xx)
            yymod = self(xx, *p)
            args += (",", xx, yymod, xx2, 'us lam?3:1:2 w l lc 3 t ""')
        if res:
            rms = np.std(y-ymod)
            gplot.mxtics().mytics().my2tics()
            # overplot residuals
            gplot.y2range('[-0.2:2]').ytics('nomirr').y2tics()
            args += (",", x, y-ymod, x2, "us lam?3:1:2 w p pt 7 ps 0.5 lc 1 axis x1y2 t 'res (%.3g \~ %.3g%%)', 0 lc 3 axis x1y2 t ''" % (rms, rms/np.mean(ymod)*100))
        gplot(*args)


def show_model(x, y, ymod, res=True):
    gplot(x, y, ymod, 'w lp pt 7 ps 0.5 t "S_i",',
          '"" us 1:3 w lp pt 6 ps 0.5 lc 3 t "S(i)"')
    if res:
        rms = np.std(y-ymod)
        gplot.mxtics().mytics().my2tics()
        # overplot residuals
        gplot.y2range('[-0.2:2]').ytics('nomirr').y2tics()
        gplot+(x, y-ymod, "w p pt 7 ps 0.5 lc 1 axis x1y2 t 'res %.3g', 0 lc 3 axis x1y2" % rms)
