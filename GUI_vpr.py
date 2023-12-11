#! /usr/bin/env python3
# -*- coding: iso-8859-1 -*-
#
## viper - Velocity and IP Estimator
## Copyright (C) Mathias Zechmeister and Jana Koehler
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along
## with this program; if not, write to the Free Software Foundation, Inc.,
## 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# GUI to run vpr.py, showing most important options

import os
from tkinter import *
from tkinter.filedialog import askopenfilename, askdirectory
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText
from utils.hbox import Help_Box

import numpy as np
import re
import sys
import vpr
from vpr import VPR
vpr.pause = print   # mainloop of the gui will pause

viperdir = os.path.dirname(os.path.realpath(__file__)) + os.sep

font_type = 'Arial'
font_size = 12
bg_frame = '#f1f1f1'    # bg color big frame
bg_color = '#e6e1e1'    # bg color small frames
bg_button = "#fdfdfd"   # bg buttons

win_width = 880     # width of GUI window
win_high = 830      # height of GUI window
fr_high = win_high-115
xy0 = 10


def set_oset(cbv, value):
    # set all values of a checkbutton-array to one value
    for v in cbv:
        v.set(value)
 
def CheckBoxes(frm, para, cbi, cbv, cols):
    if cbi != []:
        # delete old checkboxes
        for i in cbi: i.destroy()
        cbi.clear()
        cbv.clear()

    for i, o in enumerate(para):
        cbv.append(IntVar())
        yi, xi =  divmod(i, cols)
        c = ttk.Checkbutton(frm, text="  "+str(o),  variable=cbv[i])
        c.grid(row=yi+1, column=xi, sticky="nw", padx=8, pady=2)
        cbi.append(c)
        frm.grid_columnconfigure(xi, weight=1)

    return cbi, cbv  

def text_from_file(text):
     # read in help text from vpr.py argopt
     with open(viperdir+"vpr.py") as search:
         for line in search:
             line = line.rstrip()  
             if text in line:
                  infotext = (line.split('help')[1]).split("'")[1]
     return infotext


class TabFrame(Frame):
    """
    Generate Frame for the individual Tabs
    """
    def __init__(self, master, fr_text):

        self.master = master

        frm = Frame(self.master, bg=bg_frame, bd=2, relief='groove')
        frm.grid(row=0, column = 0, sticky="news",padx=10,pady=(10,10), columnspan=2)
        frm.grid_propagate(False)
        frm.columnconfigure(0, weight=1, uniform='x')
        frm.columnconfigure(1, weight=1, uniform='x')

        lb = ttk.Label(frm, text=fr_text, font=(font_type, font_size, 'bold'))
        lb.grid(row=0, column=0, sticky="nw", padx=xy0, pady=(xy0,5), columnspan = 2)

        self.frm = frm


class CBFrame(Frame):
    """
    Generate LabelFrame, including 'select all' and 'select none' Buttons
    """
    def __init__(self, master, fr_text, cbi):

        self.master = master
        self.cbi = cbi
        lfrm = LabelFrame(master, text=fr_text, bg=bg_frame, bd=2)

        bt_all = ttk.Button(lfrm, text='select all', style = 'small.TButton', command = lambda: set_oset(self.cbi,1))
        bt_all.grid(row=0, column=0, sticky="ne", padx=10, pady=(0,5), columnspan=6)
        bt_none = ttk.Button(lfrm, text='select none', style = 'small.TButton', command = lambda: set_oset(self.cbi,0))
        bt_none.grid(row=0, column=0, sticky="ne", padx=100, pady=(0,5), columnspan=6)

        self.lfrm = lfrm


class RV_plot(ttk.Frame):
    """
    Frame for RV plotting options
    """
    def __init__(self, master, e_run):

        self.filename1 = StringVar()
        self.filename2 = StringVar()
        self.cb_oset1, self.cb_oset2 = [], []		# oset checkboxes
        self.cbv_oset1, self.cbv_oset2 = [], []		# Var() of  oset checkboxes 
        self.o_rvo1, self.o_rvo2 = [], []		# order names
        self.plt_opt = '-plot rv'
        self.plt_cmp = False
        self.e_run = e_run

        self.frm_rv = TabFrame(master, 'Plot RVs').frm
        self.frm_rv.grid_rowconfigure(100, weight=1)

        # Frames
        frm_rv1 = LabelFrame(self.frm_rv, text="Input", bg=bg_frame, bd=2)
        frm_rv1.grid(row=1, column=0, sticky="news", padx=10, pady=5, ipady=5, columnspan=2)
        frm_rv1.grid_columnconfigure(1, weight=1)

        self.lfr_oset1 = CBFrame(self.frm_rv, 'oset rvo 1', self.cbv_oset1).lfrm
        self.lfr_oset1.grid(row=2, column=0, sticky="news", padx=(10,0), pady=5, ipady=5)

        self.lfr_oset2 = CBFrame(self.frm_rv, 'oset rvo 2', self.cbv_oset2).lfrm
        self.lfr_oset2.grid(row=2, column=1, sticky="news", padx=10, pady=5, ipady=5)

        lfr_cen = LabelFrame(self.frm_rv, text="center RV values", bg=bg_frame, bd=2)
        lfr_cen.grid(row=3, column=0, sticky="news", padx=(10,0), pady=5, ipady=5)

        lfr_plot = LabelFrame(self.frm_rv, text="Plotting Options", bg=bg_frame, bd=2)
        lfr_plot.grid(row=3, column=1, sticky="news", padx=10, pady=5, ipady=5)

        lfr_other = LabelFrame(self.frm_rv, text="Other", bg=bg_frame, bd=2)
        lfr_other.grid(row=4, column=0, sticky="news", padx=10, pady=5, ipady=5, columnspan=2)
        lfr_other.grid_columnconfigure(1, weight=1)

        for cc in range(0,5,1):
            self.lfr_oset1.grid_columnconfigure(cc, weight=1)
            self.lfr_oset2.grid_columnconfigure(cc, weight=1)

        # Label
        l_rv1 = ttk.Label(frm_rv1, text='rvo file 1:', font=(font_type, font_size-1, 'bold'))
        l_rv1.grid(row=0, column=0, sticky="sw", padx=xy0, pady=5)

        l_rv2 = ttk.Label(frm_rv1, text='rvo file 2:', font=(font_type, font_size-1, 'bold'))
        l_rv2.grid(row=1, column=0, sticky="nw", padx=xy0, pady=5)

        l_sort = ttk.Label(lfr_plot, text='sort by:')
        l_sort.grid(row=0, column=0, sticky="nw", padx=xy0, pady=(5,0))
        Help_Box(widget = l_sort, text = text_from_file("'-sort'"))

        l_off = ttk.Label(lfr_plot, text='offset rvo:')
        l_off.grid(row=1, column=0, sticky="nw", padx=xy0, pady=(5,0))
        Help_Box(widget = l_off, text = text_from_file("'-offset'"))

        l_aver = ttk.Label(lfr_other, text='average:')
        l_aver.grid(row=0, column=0, sticky="nw", padx=xy0, pady=(5,0))
        Help_Box(widget = l_aver, text = text_from_file("'-avg'"))

        l_out = ttk.Label(lfr_other, text='output:')
        l_out.grid(row=1, column=0, sticky="nw", padx=xy0, pady=(5,0))
        Help_Box(widget = l_out, text = text_from_file("'-save'"))

        # Entry
        self.e_rvo1 = Entry(frm_rv1, textvariable=self.filename1, width=100)
        if os.path.isfile('tmp.rvo.dat'):
            self.e_rvo1.insert(0, 'tmp.rvo.dat')
        elif os.path.isfile('tmp_rvo_par.fits'):
            self.e_rvo1.insert(0, 'tmp_rvo_par.fits')
        self.e_rvo1.bind("<Return>", (lambda event: self.update(refresh=1)))
        self.e_rvo1.grid(row=0, column=1, sticky="sew", padx=xy0, pady=(5))

        self.e_rvo2 = Entry(frm_rv1, textvariable=self.filename2, width=100)
        self.e_rvo2.insert(0, '')
        self.e_rvo2.bind("<Return>", (lambda event: self.update(refresh=2)))
        self.e_rvo2.grid(row=1, column=1, sticky="new", padx=xy0, pady=(5))

        self.e_sort = Entry(lfr_plot, width=10)
        self.e_sort.insert(0, 'BJD')
        self.e_sort.grid(row=0, column=1, sticky="nw", padx=xy0, pady=(5,0))

        self.e_offset = Entry(lfr_plot, width=10)
        self.e_offset.insert(0, 400)
        self.e_offset.bind("<Return>", (lambda event: self.update()))
        self.e_offset.grid(row=1, column=1, sticky="nw", padx=(xy0), pady=(5,0))

        self.e_out = Entry(lfr_other, width=350)
        self.e_out.insert(0, 'tmp.dat')
        self.e_out.grid(row=1, column=1, sticky="nw", padx=xy0, pady=(5,0))

        # Combobox
        self.combo_avg = ttk.Combobox(lfr_other, values=['mean', 'wmean'], width=8)
        self.combo_avg.set('wmean')
        self.combo_avg.bind('<<ComboboxSelected>>', (lambda event: self.update()))
        self.combo_avg.grid(row=0, column=1, sticky="nw", padx=xy0, pady=(5,0))
         
        # Checkbutton
        self.cb_cen = IntVar()
        l_cen = ttk.Checkbutton(lfr_cen, text="     cen RVs to zero median", variable=self.cb_cen)
        l_cen.grid(row=0, column=0, sticky="nw", padx=xy0, pady=(5,0))
        Help_Box(widget = l_cen, text = text_from_file("'-cen'"))

        self.cb_ocen = IntVar()
        l_ocen1 = ttk.Checkbutton(lfr_cen, text="     ocen rvo 1", variable=self.cb_ocen)
        l_ocen1.grid(row=1, column=0, sticky="nw", padx=xy0, pady=(5,0))
        Help_Box(widget = l_ocen1, text = text_from_file("'-ocen'"))

        self.cb_cmpocen = IntVar()
        l_ocen2 = ttk.Checkbutton(lfr_cen, text="     ocen rvo 2", variable=self.cb_cmpocen)
        l_ocen2.grid(row=2, column=0, sticky="nw", padx=xy0, pady=(5,0))
        Help_Box(widget = l_ocen2, text = text_from_file("'-cmpocen'"))

        self.cb_cen.trace("w", self.update)
        self.cb_ocen.trace("w", self.update)
        self.cb_cmpocen.trace("w", self.update)

        # Buttons
        b_rvo1 = ttk.Button(frm_rv1, text='Search data file', style='nbold.TButton', command = self.bt_rvo1)
        b_rvo1.grid(row=0, column=2, sticky="ne", padx=xy0, pady=(0,3))

        b_rvo2 = ttk.Button(frm_rv1,text='Search data file', style='nbold.TButton', command = self.bt_rvo2)
        b_rvo2.grid(row=1, column=2, sticky="ne", padx=xy0, pady=(3))

        b_swap = ttk.Button(frm_rv1, text='Swap', style='nbold.TButton', command = self.bt_swap)
        b_swap.grid(row=2, column=1, sticky="ne", padx=xy0, pady = 3, rowspan=2)

        b_cmp = ttk.Button(frm_rv1, text='Compare', style='bold.TButton', command = lambda: self.show_plot(cmp=True))
        b_cmp.grid(row=2, column=2, sticky="ne", padx=xy0, pady = 3)

        b_save = ttk.Button(lfr_other, text='Save', style='nbold.TButton', command = lambda: self.show_plot('-save'))
        b_save.grid(row=1, column=2, sticky="ne", padx=xy0)

        b_nrvo = ttk.Button(self.frm_rv, text='Plot no-rv', style='bold.TButton', command = lambda: self.show_plot('-plot nrvo'))
        b_nrvo.grid(row=100, column=0, sticky="se", padx=20, pady=(5,20), columnspan=10)

        b_rvo = ttk.Button(self.frm_rv, text='Plot o-rv', style='bold.TButton', command = lambda: self.show_plot('-plot rvo'))
        b_rvo.grid(row=100, column=0, sticky="se", padx=180, pady=(5,20), columnspan=10)

        b_rvbjd = ttk.Button(self.frm_rv, text='Plot BJD-RV', style='bold.TButton', command = lambda: self.show_plot())
        b_rvbjd.grid(row=100, column=0, sticky="se", padx=340, pady=(5,20), columnspan=10)

        if sys.argv[1:]:
            self.e_rvo1.delete(0, END)
            self.e_rvo1.insert(0, sys.argv[1])

        #if os.path.isfile(self.e_rvo1.get()):
        self.update(refresh=1, plot=0, destroy = 0)

    def bt_rvo1(self):
        file1 = askopenfilename()
        if file1:
            self.filename1.set(file1)
            self.update(refresh=1)

    def bt_rvo2(self):
        file1 = askopenfilename()
        if file1:
            self.filename2.set(file1)
            self.update(refresh=2)

    def bt_swap(self):
        f1 = self.filename1.get()
        f2 = self.filename2.get()
        self.filename1.set(f2)
        self.filename2.set(f1)
        self.update(refresh=12)

    def get_orders(self, file):
        # get all available orders from RV file

        Afull = VPR(tag=file).Afull 
        colnames = Afull.dtype.names

        onames_all = np.array([col for col in colnames if col.startswith('rv')], dtype='O')
        orders_all, chunks_all = [*np.array([[*map(int, oname[2:].split('-'))] for oname in onames_all]).T, None][0:2]
        return orders_all

    def refresh_oset(self, num):
        # refresh oset checkboxes (when using new RV file)
        if '1' in num:
            if self.e_rvo1.get():
                self.o_rvo1 = self.get_orders(self.e_rvo1.get())
                self.cb_oset1, self.cbv_oset1 = CheckBoxes(self.lfr_oset1, self.o_rvo1, self.cb_oset1, self.cbv_oset1, 6)
                set_oset(self.cbv_oset1, 1)
                for c in self.cbv_oset1: c.trace("w", self.update)
        if '2' in num:
            if self.e_rvo2.get():
                self.o_rvo2 = self.get_orders(self.e_rvo2.get())
                self.cb_oset2, self.cbv_oset2 = CheckBoxes(self.lfr_oset2, self.o_rvo2, self.cb_oset2, self.cbv_oset2, 6)
                set_oset(self.cbv_oset2, 1)
                for c in self.cbv_oset2: c.trace("w", self.update)

    def update(self, *args, refresh = 0, plot = 1, destroy = 1):
        try:
            self.refresh_oset(str(refresh))
            if plot: self.show_plot(self.plt_opt, cmp=self.plt_cmp)
            if destroy and hasattr(self, 'emsg'): self.emsg.destroy()
        except:
            self.emsg = ttk.Label(self.frm_rv, text="Problem using input file. Check if selected file exists and is not empty.", foreground='red')
            self.emsg.grid(row=0, column=0, sticky="nw", padx=(200,0), pady=(xy0,5), columnspan = 10)

    def show_plot(self, args='-plot rv', cmp=False):
        # plot selected options
        self.plt_cmp = cmp
        set_oset1 = np.array([o.get() for o in self.cbv_oset1])
        set_oset2 = np.array([o.get() for o in self.cbv_oset2])

        str_arg = self.e_rvo1.get() + ' -avg ' + self.combo_avg.get()

        if self.cb_cen.get(): str_arg += " -cen "
        if self.cb_ocen.get(): str_arg += " -ocen "
        if self.e_sort.get(): str_arg += " -sort " + str(self.e_sort.get())
        if self.e_offset.get(): str_arg += " -offset " + str(self.e_offset.get())

        if set_oset1.any():
            str_arg += " -oset " 
            for o in self.o_rvo1[set_oset1==1]: str_arg += str(o) + ","

        if cmp:
            if self.e_rvo2.get(): str_arg += " -cmp " + self.e_rvo2.get()
            if self.cb_cmpocen.get(): str_arg += " -cmpocen"
            if set_oset2.any():
                str_arg += " -cmposet " 
                for o in self.o_rvo2[set_oset2==1]: str_arg += str(o) + ","

        if '-save' in str(args):    
            str_arg += ' -save ' + self.e_out.get()
            args= '-plot rv'
        
        if args:
            self.plt_opt = args
            str_arg += " " + args
 
        if self.cbv_oset1 != []:
            self.e_run.delete('0.0', END)
            self.e_run.insert(INSERT,"python3 vpr.py "+str_arg)
            vpr.run(str_arg.split())


class Parameter(ttk.Frame):
    """
    Frame for Parameter plotting options
    """
    def __init__(self, master, e_run):

        self.cb_parx, self.cb_pary = StringVar(), StringVar()

        self.frm_par = TabFrame(master, 'Plot parameter').frm
        self.frm_par.grid_rowconfigure(100, weight=1)
        self.e_run = e_run

        # Frames
        lpar = LabelFrame(self.frm_par, text='Input', bg=bg_frame, bd=2)
        lpar.grid(row=1, column=0, sticky="news", padx=10, pady=5, ipady=5, columnspan=10)
        lpar.grid_columnconfigure(1, weight=1)

        # Labels
        l_parf = ttk.Label(lpar, text='parfile:', font=(font_type, font_size-1, 'bold'))
        l_parf.grid(row=1, column=0, sticky="nw", padx=xy0, pady=5)

        l_x = ttk.Label(self.frm_par, text='Parameters x-axis', font=(font_type, font_size-1, 'bold'))
        l_x.grid(row=3, column=0, sticky="nw", padx=xy0, pady=xy0, columnspan = 5)

        l_y = ttk.Label(self.frm_par, text='Parameters y-axis', font=(font_type, font_size-1, 'bold'))
        l_y.grid(row=5, column=0, sticky="nw", padx=xy0, pady=xy0, columnspan = 5)

        # Entry
        self.e_parf = Entry(lpar, width=120)
        if os.path.isfile('tmp.par.dat'):
            self.e_parf.insert(0, 'tmp.par.dat')
        elif os.path.isfile('tmp_rvo_par.fits'):
            self.e_parf.insert(0, 'tmp_rvo_par.fits')
        self.e_parf.bind("<Return>", (lambda event: self.update()))
        self.e_parf.grid(row=1, column=1, sticky="new", padx=xy0, pady=(5))

        # Buttons
        b_fpar = ttk.Button(lpar,text='Search data file', style='nbold.TButton', command = lambda: self.bt_file())
        b_fpar.grid(row=1, column=2, sticky="ne", padx=10, pady=0)

        b_par = ttk.Button(self.frm_par, text='Plot par', style='bold.TButton', command = lambda: self.show_plot())
        b_par.grid(row=100, column=0, sticky="se", padx=20, pady=20, columnspan=10)

      #  if os.path.isfile(self.e_parf.get()):
        self.update(plot=0, destroy = 0)

    def bt_file(self):
        edir = askopenfilename()
        if edir:
            self.e_parf.delete(0, END)
            self.e_parf.insert(0, edir)
            self.update()

    def update(self, plot=1, destroy = 1):
        try:
            self.get_parameters()
            if plot: self.show_plot()
            if destroy and hasattr(self, 'emsg'): self.emsg.destroy()
        except:
            self.emsg = ttk.Label(self.frm_par, text="Problem using input file. Check if selected file exists and is not empty.", foreground='red')
            self.emsg.grid(row=0, column=0, sticky="nw", padx=(200,0), pady=(xy0,5), columnspan = 10)

    def get_parameters(self):
        parfile = self.e_parf.get()#.split('.')[0]

        if parfile.endswith('fits'):
            from astropy.io import fits
            hdu = fits.open(parfile, dtype=None)
            par = hdu[2].data.view(np.recarray) 
        else:
            par = np.genfromtxt(parfile, dtype=None, names=True,
                            deletechars='',   # to keep the dash for chunks
                            encoding=None).view(np.recarray)

        colnames = par.dtype.names[:-1]

        # parameters x axis
        l = LabelFrame(self.frm_par, text='', bg=bg_frame, bd=2)
        l.grid(row=4, column=0, sticky="news", padx=10, pady=5, ipady=5)

        for i, c in enumerate(colnames[:3]):
            cc = ttk.Radiobutton(l, text="  "+str(c),  variable=self.cb_parx, value=c)
            cc.grid(row=i, column=0, sticky="new", padx=15, pady=2)
        self.cb_parx.set('n')

        # parameters y axis
        lfr_par = []
        par_groups = ['rv', 'norm', 'wave', 'ip', 'atm', 'bkg']
        for i, cn in enumerate(par_groups):
            yi, xi =  divmod(i, 6)
            l = LabelFrame(self.frm_par, text=cn, bg=bg_frame, bd=2)        
            l.grid(row=6+yi, column=xi, sticky="news", padx=10, pady=5, ipady=5)
            self.frm_par.grid_columnconfigure(i, weight=1, uniform='x')
            lfr_par.append(l)

        for i, c in enumerate(colnames[4:][::2]):
            pos = par_groups.index(re.sub(r'[0-9]', '', c))
            cc = ttk.Radiobutton(lfr_par[pos], text="  "+str(c), variable=self.cb_pary, value=c)
            cc.grid(row=i+5, column=0, sticky="new", padx=15, pady=2)
        self.cb_pary.set('ip0')
    
        self.cb_parx.trace("w", lambda *args: self.show_plot())
        self.cb_pary.trace("w", lambda *args: self.show_plot())

    def show_plot(self):
        str_arg = self.e_parf.get().split('.par')[0] + " -plot par "  

        str_arg += ' -parcolx '+str(self.cb_parx.get())
        str_arg += ' -parcoly '+str(self.cb_pary.get())

        self.e_run.delete('0.0', END)
        self.e_run.insert(INSERT,"python3 vpr.py "+str_arg)
        vpr.run(str_arg.split())


class Residuals(ttk.Frame):
    '''
    Frame for Residual plotting options
    '''
    def __init__(self, master, e_run):

        # Arrays checkbuttons + variable for nset and oset
        self.cb_nset, self.cbv_nset = [], []
        self.cb_oset, self.cbv_oset = [], []
        self.e_run = e_run

        # Frames
        self.frm_res = TabFrame(master, 'Plot residuals').frm
        self.frm_res.grid_rowconfigure(100, weight=1)

        lres = LabelFrame(self.frm_res, text='Input', bg=bg_frame, bd=2)
        lres.grid(row=1, column=0, sticky="news", padx=10, pady=5, ipady=5, columnspan=2)
        lres.grid_columnconfigure(1, weight=1)

        self.lreso = CBFrame(self.frm_res, 'oset', self.cbv_oset).lfrm
        self.lreso.grid(row=7, column=1, sticky="news", padx=(5,10), pady=5, ipady=5)

        self.lresn = CBFrame(self.frm_res, 'nset', self.cbv_nset).lfrm
        self.lresn.grid(row=7, column=0, sticky="news", padx=(10,5), pady=5, ipady=5)

        # Labels
        l_info = ttk.Label(self.frm_res, text='For large number of nset, it is recommended to first chose one single order before using "select all" on nset.', foreground='red4')
        l_info.grid(row=2, column=0, sticky="nw", padx=xy0, pady=xy0, columnspan=2)

        l_dir = ttk.Label(lres, text='directory:', font=(font_type, font_size-1, 'bold'))
        l_dir.grid(row=1, column=0, sticky="nw", padx=xy0, pady=3)

        l_ressep = ttk.Label(lres, text='offset res:', background=bg_frame)
        l_ressep.grid(row=2, column=0, sticky="nw", padx=xy0, pady=3)
        Help_Box(widget = l_ressep, text = text_from_file("'-ressep'"))

        # Entry
        self.e_dir = Entry(lres, width = 60)
        self.e_dir.insert(0, 'res')
        self.e_dir.bind("<Return>", (lambda event: self.update()))
        self.e_dir.grid(row=1, column=1, sticky="new", padx=xy0, pady=3)

        self.e_ressep = Entry(lres, width = 20)
        self.e_ressep.insert(0, 5)
        self.e_ressep.bind("<Return>", (lambda event: self.show_plot()))
        self.e_ressep.grid(row=2, column=1, sticky="nw", padx=xy0, pady=3)

        # Buttons
        b_fres = ttk.Button(lres, text='Search directory', style='nbold.TButton', command = lambda: self.bt_file())
        b_fres.grid(row=1, column=2, sticky="ne", padx=10, pady=0)

        b_res = ttk.Button(self.frm_res, text='Plot res', style='bold.TButton', command = lambda: self.show_plot())
        b_res.grid(row=100, column=0, sticky="se", padx=20, pady=20, columnspan=10)

        #if os.path.isdir(self.e_dir.get()):        
        self.update(plot=0, destroy = 0)

    def bt_file(self):
        edir = askdirectory()
        if edir:
            self.e_dir.delete(0, END)
            self.e_dir.insert(0, edir)
            self.update()

    def sel_res(self):
        list_dat = np.array(sorted(os.listdir(self.e_dir.get())))
        num, orders = [], []
        for dat in list_dat:
            no = (dat.split(".")[0]).split("_")
            num.append(int(no[0]))
            orders.append(int(no[1]))

        self.num = np.unique(num)
        self.orders = np.unique(orders)

        self.cb_nset, self.cbv_nset = CheckBoxes(self.lresn, self.num, self.cb_nset, self.cbv_nset, 5)
        self.cb_oset, self.cbv_oset = CheckBoxes(self.lreso, self.orders, self.cb_oset, self.cbv_oset, 5)

        set_oset(self.cbv_oset, 1)
        for c in self.cbv_oset: c.trace("w", lambda *args: self.show_plot())

        self.cbv_nset[0].set(1)
        for c in self.cbv_nset: c.trace("w", lambda *args: self.show_plot())   

    def update(self, plot=1, destroy = 1):
        try:
            self.sel_res()
            if plot: self.show_plot()
            if destroy and hasattr(self, 'emsg'): self.emsg.destroy()
        except:
            self.emsg = ttk.Label(self.frm_res, text="Problem using input directory. Check if selected directory exists and is not empty.", foreground='red')
            self.emsg.grid(row=0, column=0, sticky="nw", padx=(200,0), pady=(xy0,5), columnspan = 10)

    def show_plot(self):
        str_arg = " -res "+str(self.e_dir.get())

        set_nset = np.array([o.get() for o in self.cbv_nset])
        set_oset = np.array([o.get() for o in self.cbv_oset])

        if set_nset.any():
            str_arg += " -nset " + str((self.num[set_nset==1]).tolist()).replace(' ','') 

        if set_oset.any():
            str_arg += " -oset " + str((self.orders[set_oset==1]).tolist()).replace(' ','')  

        if self.e_ressep.get(): str_arg += " -ressep "+ self.e_ressep.get()

        self.e_run.delete('0.0', END)
        self.e_run.insert(INSERT,"python3 vpr.py "+str_arg)
        vpr.run(str_arg.split())


def main():
        win = Tk()
        win.title('Gui VPR')
        win.geometry("{}x{}".format(win_width, win_high))
        win.configure(background=ttk.Style().lookup('TFrame', 'background'), padx=3, pady=6)  

        win.grid_propagate(False)
        win.columnconfigure(0, weight=1)
        win.rowconfigure(0, weight=10, uniform=1, minsize=fr_high)
        win.rowconfigure(1, weight=2, uniform=1)
        win.rowconfigure(2, weight=10, uniform=1)

        s = ttk.Style()
        s.configure("TCheckbutton", background=bg_frame, bd=0, highlightthickness=0)     
        s.configure("bold.TButton", padding=2, background=bg_button, font=(font_type,font_size,'bold'), borderwidth=2, width=15)
        s.configure("nbold.TButton", padding=2, background=bg_button, font=(font_type,font_size,''), borderwidth=2, width=15)
        s.configure("small.TButton", padding=1, background=bg_button, font=(font_type,font_size-1,''), borderwidth=1)       
        s.configure('TRadiobutton', background=bg_frame)
        s.configure('TFrame', background=bg_frame)
        s.configure('TLabel', background=bg_frame)
        s.configure('lefttab.TNotebook', tabposition='nw')
        s.configure('TNotebook.Tab', font=(font_type, font_size+1, 'bold'))
        s.map('TNotebook.Tab', background= [("selected", bg_color)])

        tabControl = ttk.Notebook(win, style='lefttab.TNotebook')
        tab_rv = ttk.Frame(tabControl, style='TFrame')
        tab_par = ttk.Frame(tabControl, style='TFrame')
        tab_res = ttk.Frame(tabControl, style='TFrame')

        for tab in ([tab_rv, tab_par, tab_res]):
            tab.grid_propagate(False)
            tab.columnconfigure(0, weight=1)
            tab.rowconfigure(0, weight=1)

        tabControl.add(tab_rv, text ='  RV plots  ')
        tabControl.add(tab_par, text ='  Parameters  ')
        tabControl.add(tab_res, text ='  Residuals  ')
        tabControl.grid(row=0, column = 0, sticky="news",padx=3,pady=(3,5), columnspan=2)

        lrun = Label(win, text='Current command:')
        lrun.grid(row=1, column = 0, sticky="w",padx=10,pady=0)
        e_run = ScrolledText(win, background='#f0f0f0', width=200)
        e_run.grid(row=2, column = 0, sticky="nsew",padx=(10,50),pady=(5,10))

        # set up the Frames for all Tabs
        RV_plot(tab_rv, e_run)
        Parameter(tab_par, e_run)
        Residuals(tab_res, e_run)

        b_exit = ttk.Button(master=win, text='EXIT', style='bold.TButton', width=10, command = lambda: exit())
        b_exit.grid(row=2, column = 1, sticky="se",padx=10,pady=(5,10))

        win.mainloop()

if __name__ == '__main__':
    main()
