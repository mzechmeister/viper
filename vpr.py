#! /usr/bin/env python3

# ./vpr.py

import numpy as np

from pause import pause
from gplot import gplot
gplot.tmp = '$'

def plot_rvo(rv=None, e_rv=None, file=None):
    #gplot('"%s" matrix every 2::1 us (($1-1)/2+18):3 w p pt 7' % file)
    if file:
        mat = np.genfromtxt(file)
        bjd, RV, e_RV = mat[:3]
        rv, e_rv = mat[3::2], mat[4::2]

    ii = np.isfinite(e_rv)
    RV = np.mean(rv[ii]) 
    e_RV = np.std(rv[ii])/(ii.sum()-1)**0.5

    gplot.xlabel('"Order index"')
    gplot.ylabel('"RV [km/s]"')
    gplot(rv, e_rv, 'us 0:1:2 w e pt 7 t "%s", RV=%s, e_RV=%s, RV lc 3 t "RV = %.5f +/- %.5f km/s", RV+e_RV lc 3 dt 2 t "", RV-e_RV lc 3 dt 2 t ""' % (file, RV,e_RV,RV,e_RV))
    pause()

def plot_rv(rv=None, e_rv=None, file=None):
    gplot('"%s" us 1:2:3 w e pt 7' % file)
    pause()

if __name__ == "__main__":
    plot_rv(file="tmp.rvo.dat")
    plot_rvo(file="tmp.rvo.dat")

