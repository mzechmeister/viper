#! /usr/bin/env python3

# ./vpr.py

import numpy as np

from pause import pause
from gplot import gplot
gplot.tmp = '$'

def plot_RV(rv=None, e_rv=None, file=None):
    gplot('"%s" us 1:2:3 w e pt 7' % file)
    pause()

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


class VPR():
    def __init__(self, file):
        self.A = np.genfromtxt(file, usecols=(0, 1, 2), dtype=[('bjd','float'), ('RV',float), ('e_RV',float)]).view(np.recarray)
        mat = np.genfromtxt(file)
        self.rv, self.e_rv = mat.T[5::2], mat.T[6::2]
    def plot_rv(self, o=None, n=None):
        A = self.A
        gplot.var(n=1, N=len(self.rv.T))
        gplot.bind('")" "n = n>=N? N : n+1; repl"')
        gplot.bind('"(" "n = n<=1? 1 : n-1; repl"')
        gplot.xlabel("'order o'")
        gplot.ylabel("'RV_{n,o} -- RV_{n}  [m/s]'")
        gplot('for [n=1:N]', (self.rv-A.RV).T*1000, self.e_rv.T*1000, 'us ($0-0.25+0.5*n/N):(column(n)):(column(n+N)) w e pt 6 lc "light-grey" t "", "" us ($0-0.25+0.5*n/N):(column(n)):(column(n+N)) w e pt 7 lc 1 t "n = ".n')
        pause()
        
if __name__ == "__main__":
    VPR("tmp.rvo.dat").plot_rv(o=1)
    plot_RV(file="tmp.rvo.dat")
    plot_rvo(file="tmp.rvo.dat")

