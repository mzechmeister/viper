#! /usr/bin/env python3

# ./vpr.py

import argparse

import numpy as np

from pause import pause
from gplot import gplot

def plot_RV(file):
    gplot.mxtics().mytics()
    gplot.xlabel("'BJD - 2 450 000'")
    gplot.ylabel("'RV [m/s]'")
    gplot.key("title '%s'" % file)
    gplot('"%s" us ($1-2450000):"RV":(sprintf("%%s\\nn: %%d\\nBJD: %%.6f\\nRV: %%f +/- %%f", strcol("filename"),$0+1,$1,$2,$3)) with labels  hypertext point pt 0 t"", "" us ($1-2450000):"RV":"e_RV" w e pt 7 lc 7' % file)
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
    gplot.ylabel('"RV [m/s]"')
    gplot(rv, e_rv, 'us 0:1:2 w e pt 7 t "%s", RV=%s, e_RV=%s, RV lc 3 t "RV = %.5f +/- %.5f m/s", RV+e_RV lc 3 dt 2 t "", RV-e_RV lc 3 dt 2 t ""' % (file, RV,e_RV,RV,e_RV))
    pause()


class VPR():
    def __init__(self, tag):
        self.tag = tag
        self.file = file = tag + '.rvo.dat'
        self.A = np.genfromtxt(file, usecols=range(-1,4), dtype=None, names=True, encoding=None).view(np.recarray)
        mat = np.genfromtxt(file, skip_header=1)
        self.rv, self.e_rv = mat.T[4:-1:2], mat.T[5::2]
        orders = np.genfromtxt(self.file, names=True).dtype.names[4:-1:2]
        self.orders = [int(o[2:]) for o in orders]

    def plot_RV(self):
        plot_RV(self.file)

    def plot_rv(self, o=None, n=None):
        A = self.A
        gplot.var(n=1, N=len(self.rv.T))
        gplot.key('title "%s"'%self.tag)
        gplot.bind('")" "n = n>=N? N : n+1; repl"')
        gplot.bind('"(" "n = n<=1? 1 : n-1; repl"')
        gplot.xlabel("'order o'")
        gplot.ylabel("'RV_{n,o} -- RV_{n}  [m/s]'")
        gplot.mxtics().mytics()

        gplot('for [n=1:N]', self.orders, (self.rv-A.RV).T, self.e_rv.T, 'us ($1-0.25+0.5*n/N):(column(1+n)):(column(1+n+N)) w e pt 6 lc "light-grey" t "", "" us ($1-0.25+0.5*n/N):(column(1+n)):(column(1+n+N)) w e pt 6 lc 1 t "RV_{".n.",o} -- RV_{".n."}",', A.BJD, A.RV+400, A.e_RV, A.filename, ' us 1:2:(sprintf("%s\\nn: %d\\nBJD: %.6f\\nRV: %f +/- %f",strcol(4),$0+1,$1,$2,$3)) w labels hypertext point pt 0 axis x2y1 t "", "" us 1:2:3 w e lc 7 pt 7 axis x2y1 t "RV_n", "" us 1:($2/($0+1==n)):3 w e lc 1 pt 7 axis x2y1 t "RV_{".n."}"')
        pause()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Analyse viper RVs', add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    argopt = parser.add_argument   # function short cut
    argopt('tag', nargs='?', help='tag', default='tmp', type=str)

    args = parser.parse_args()

    vpr = VPR(args.tag)
    vpr.plot_RV()
    vpr.plot_rv(o=1)
