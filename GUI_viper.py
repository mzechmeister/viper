#! /usr/bin/env python3
# -*- coding: iso-8859-1 -*-
# Licensed under a GPLv3 style license - see LICENSE

# GUI to start viper.py, showing most important options

import os
import sys
import glob
from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter import ttk
import numpy as np
import importlib
import configparser
from tkinter.scrolledtext import ScrolledText
from utils.hbox import Help_Box
from model import IPs

viperdir = os.path.dirname(os.path.realpath(__file__)) + os.sep


font_type = 'Arial'
font_size = 12
bg_frame = '#f1f1f1'    # bg color big frame
bg_color = '#e6e1e1'    # bg color small frames
bg_button = "#fdfdfd"   # bg buttons

win_width = 1060     # width of GUI window	
win_high = 800      # height of GUI window	
xy0 = 20
y1 = 3
x1 = 10


def text_from_file(text):
     with open(viperdir+"viper.py") as search:
         for line in search:
             line = line.rstrip()  
             if text in line:
                  infotext = (line.split('help')[1]).split("'")[1]
     return infotext


class GUI_viper():
    """
    Build up Frames and functions
    """
    def __init__(self, master, configs):

        self.cb_atm = []       # Checkboxes telluric molecules 
        self.cbv_atm = []      # Variables telluric molecules
        self.cb_lookctpl = IntVar()

        self.master = master
        self.configs = configs

        self.fr_input()
        self.fr_para()
        self.fr_plot()

        self.Update_inst()
        self.Update_tell()
        self.Update_ctpl()	

        lrun = ttk.Label(master, text='Current command:', background=bg_color)
        lrun.grid(row=2, column = 1, sticky="nw", padx=0, pady=(10,10), columnspan=6)
        self.e_run = ScrolledText(master, background=bg_frame)
        self.e_run.grid(row=3, column = 1, sticky="news", padx=(0, xy0), pady=(0,10))

        b_exit = ttk.Button(master, text='EXIT', style='bold.TButton', command = lambda: exit())
        b_exit.grid(row=4, column=1, sticky="se", padx=(0,xy0), pady = (0,20))

        b_go = ttk.Button(master, text='Start', style='bold.TButton', command = self.bt_start)
        b_go.grid(row=4, column=1, sticky="se", padx=(0,180), pady = (0,20))
 
    def fr_input(self):
        # Frame for data input
        fr1 = Frame(self.master, bg=bg_frame, bd=2, relief='groove')
        fr1.grid(row=0, column = 0, sticky="news", padx=20, pady=(xy0,6), ipady=10, columnspan=2)
        fr1.grid_propagate(False)
        fr1.grid_columnconfigure(2, weight=1)
        fr1.grid_columnconfigure(4, weight=1)
        fr1.grid_rowconfigure(0, weight=1)
        fr1.grid_rowconfigure(5, weight=1)

        # Label
        l_data = ttk.Label(fr1, text='data files', font=(font_type, font_size, 'bold'))
        l_data.grid(row=1, column=0, sticky="nw", padx=x1, pady=y1)

        l_tpl = ttk.Label(fr1, text='template file', font=(font_type, font_size, 'bold'))
        l_tpl.grid(row=2, column=0, sticky="nw", padx=x1, pady=y1)

        l_inst = ttk.Label(fr1, text='spectrograph:')
        l_inst.grid(row=5, column=0, sticky="nw", padx=x1, pady=y1)
        Help_Box(widget = l_inst, text = text_from_file("'-inst'"))

        l_targ = ttk.Label(fr1, text='targ:')
        l_targ.grid(row=5, column=3, sticky="nw", padx=xy0, pady=y1)
        Help_Box(widget = l_targ, text = text_from_file("'-targ'"))

        # Entry
        self.e_dat = Entry(fr1, width = 200)
        self.e_dat.insert(0, 'data/TLS/hd189733/*')
        self.e_dat.grid(row=1, column=1, sticky="nw", padx=x1, pady=y1, columnspan=6)

        self.e_tpl = Entry(fr1, width = 200)
        self.e_tpl.insert(0, 'data/TLS/Deconv/HARPS*fits')
        self.e_tpl.grid(row=2, column=1, sticky="nw", padx=x1, pady=y1, columnspan=6)

        self.e_cell = Entry(fr1, width = 200)
        self.e_cell.grid(row=3, column=1, sticky="nw", padx=x1, pady=y1, columnspan=6)

        self.e_flag = Entry(fr1, width = 200)
        self.e_flag.insert(0, '')
        self.e_flag.grid(row=4, column=1, sticky="nw", padx=x1, pady=y1, columnspan=6)

        self.e_targ = Entry(fr1, width = 200)
        self.e_targ.insert(0, self.configs.get('targ',''))
        self.e_targ.grid(row=5, column=4, sticky="nw", padx=x1, pady=y1)

        # Checkboxes
        self.cb_cell = IntVar()
        self.cb_flagfile = IntVar()

        l_cell = ttk.Checkbutton(fr1, text="   Cell file:", variable=self.cb_cell)
        l_cell.grid(row=3, column=0, sticky="nw", padx=x1, pady=y1)
        Help_Box(widget = l_cell, text = text_from_file("'-fts'") + " If checkbox is unset, no use of cell FTS for the modelling.")
        self.cb_cell.set(not self.configs.get('nocell', 0)) 

        l_flag = ttk.Checkbutton(fr1, text="   flag file:", variable=self.cb_flagfile)
        l_flag.grid(row=4, column=0, sticky="nw", padx=x1, pady=y1)
        Help_Box(widget = l_flag, text = text_from_file("'-flagfile'"))

        insts = sorted([os.path.basename(i)[5:-3] for i in glob.glob(viperdir+'inst/inst_*.py')])
        self.combo_inst = ttk.Combobox(fr1, values=insts, width=15)
        self.combo_inst.set(self.configs.get('inst','TLS'))
        self.combo_inst.grid(row=5, column=1, sticky="nw", padx=x1, pady=y1)
        self.combo_inst.bind('<<ComboboxSelected>>', lambda event: self.Update_inst())

        # Buttons
        b_dat = ttk.Button(fr1, text='Search data file', style='nbold.TButton', command = lambda: self.bt_file(self.e_dat))
        b_dat.grid(row=1, column=7, sticky="nw", padx=xy0)

        b_tpl = ttk.Button(fr1,text='Search tpl file', style='nbold.TButton', command = lambda: self.bt_file(self.e_tpl))
        b_tpl.grid(row=2, column=7, sticky="nw", padx=xy0)

        b_cell = ttk.Button(fr1,text='Search Cell file', style='nbold.TButton', command = lambda: self.bt_file(self.e_cell))
        b_cell.grid(row=3, column=7, sticky="nw", padx=xy0)

        b_flag = ttk.Button(fr1,text='Search flag file', style='nbold.TButton', command = lambda: self.bt_file(self.e_flag))
        b_flag.grid(row=4, column=7, sticky="nw", padx=xy0)

    def fr_para(self):
        # Frame for parameter selection data reduction
        fr2 = Frame(self.master, bg=bg_frame, bd=2, relief='groove')
        fr2.grid(row=1, column = 0, sticky="news", padx=(xy0,6), pady=(0,20), rowspan=4)
        fr2.grid_propagate(False)
        fr2.grid_columnconfigure(0, weight=1, uniform='x')
        fr2.grid_columnconfigure(1, weight=1, uniform='x')

        lfr_data = LabelFrame(fr2, text="Data", bg=bg_frame, bd=2)
        lfr_data.grid(row=1, column=0, sticky="news", padx=(10,0), pady=y1, ipady=5)

        lfr_model = LabelFrame(fr2, text="Model Setup", bg=bg_frame, bd=2)
        lfr_model.grid(row=1, column=1, sticky="news", padx=(10,10), pady=y1, ipady=5)

        lfr_tpl = LabelFrame(fr2, text="Template", bg=bg_frame, bd=2)
        lfr_tpl.grid(row=2, column=1, sticky="news", padx=(10,10), pady=y1, ipady=5)

        lfr_stat = LabelFrame(fr2, text="Fit Settings", bg=bg_frame, bd=2)
        lfr_stat.grid(row=2, column=0, sticky="news", padx=(10,0), pady=y1, ipady=5)

        self.lfr_tell = LabelFrame(fr2, text="Tellurics", bg=bg_frame, bd=2)
        self.lfr_tell.grid(row=3, column=0, sticky="news", padx=(10,0), pady=y1, ipady=5, columnspan=1)

        self.lfr_ctpl = LabelFrame(fr2, text="Create Template", bg=bg_frame, bd=2)
        self.lfr_ctpl.grid(row=3, column=1, sticky="news", padx=(10,10), pady=y1, ipady=5, columnspan=1)

        lfr_out = LabelFrame(fr2, text="Output", bg=bg_frame, bd=2)
        lfr_out.grid(row=4, column=0, sticky="news", padx=(10,10), pady=y1, ipady=5, columnspan=2)

        for lfr in [lfr_data, lfr_tpl, lfr_model, lfr_stat, self.lfr_tell]:
            lfr.grid_columnconfigure(1, weight=1)

        for c in range(0,5,1):
            self.lfr_tell.grid_columnconfigure(c, weight=1)

        for c in range(1,4,1):
            lfr_out.grid_columnconfigure(c, weight=1)

        # Label
        l_opt = ttk.Label(fr2, text='Options data reduction', font=(font_type, font_size, 'bold'))
        l_opt.grid(row=0, column=0, sticky="nw", padx=(xy0,0), pady=(10,y1), columnspan=3)

        self.Label_list(lfr_data, ['nset', 'oset', 'chunks', 'vcut', 'iset'])
        self.Label_list(lfr_model, ['ip', 'iphs', 'deg_norm', 'deg_wave', 'deg_bkg'])
        self.Label_list(lfr_tpl, ['rv_guess', 'oversampling'])
        self.Label_list(lfr_stat, ['kapsig'])
        self.Label_list(lfr_out, ['tag', 'output_format'])

        self.l_tell = ttk.Label(self.lfr_tell,text='telluric:')
        self.l_tell.grid(row=0, column=0, sticky="news", padx=(xy0,5), pady=y1)
        Help_Box(widget = self.l_tell, text = text_from_file("'-telluric'"))

        # Entry
        self.e_nset = Entry(lfr_data)
        self.e_nset.insert(0, ':4')
        self.e_nset.grid(row=0, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

        self.e_oset = Entry(lfr_data)
        self.e_oset.insert(0, self.configs.get('oset','20'))
        self.e_oset.grid(row=1, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

        self.e_ch = Entry(lfr_data)
        self.e_ch.insert(0, self.configs.get('chunks','1'))
        self.e_ch.grid(row=2, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

        self.e_vcut = Entry(lfr_data)
        self.e_vcut.insert(0, self.configs.get('vcut','100'))
        self.e_vcut.grid(row=3, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

        self.e_iset = Entry(lfr_data)
        self.e_iset.grid(row=4, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0)) 

        self.e_iphs = Entry(lfr_model)
        self.e_iphs.insert(0, self.configs.get('iphs','50'))
        self.e_iphs.grid(row=1, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

        self.e_deg_norm = Entry(lfr_model)
        self.e_deg_norm.insert(0, self.configs.get('deg_norm','3'))
        self.e_deg_norm.grid(row=2, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

        self.e_deg_wave = Entry(lfr_model)
        self.e_deg_wave.insert(0, self.configs.get('deg_wave','3'))
        self.e_deg_wave.grid(row=3, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

        self.e_deg_bkg = Entry(lfr_model)
        self.e_deg_bkg.insert(0, self.configs.get('deg_bkg','1'))
        self.e_deg_bkg.grid(row=4, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

        self.e_vg = Entry(lfr_tpl)
        self.e_vg.insert(0, self.configs.get('rv_guess','1'))
        self.e_vg.grid(row=0, column=1, sticky="nw", padx=(x1,xy0), pady=(0,0))

        self.e_overs = Entry(lfr_tpl)
        self.e_overs.insert(0, self.configs.get('oversampling','1'))
        self.e_overs.grid(row=1, column=1, sticky="nw", padx=(x1,xy0), pady=0)

        self.e_kapsig = Entry(lfr_stat)
        self.e_kapsig.insert(0, self.configs.get('kapsig','4.5'))
        self.e_kapsig.grid(row=0, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

        self.e_tag = Entry(lfr_out)
        self.e_tag.insert(0, self.configs.get('tag','tmp'))
        self.e_tag.grid(row=0, column=1, sticky="news", padx=(xy0,10), pady=y1, columnspan=3)

        # Checkboxes:
        self.cb_wgt = IntVar()
        self.cb_createtpl = IntVar()
        self.cb_tellshift = IntVar()
        self.cb_format = [IntVar(), IntVar(), IntVar()]

        l_wei = ttk.Checkbutton(lfr_stat, text="     weighted error", variable=self.cb_wgt)
        l_wei.grid(row=1, column=0, sticky="nw", padx=(xy0,x1), pady=y1, columnspan=2)
        self.cb_wgt.set(self.configs.get('wgt', 0))
        Help_Box(widget = l_wei, text = text_from_file("'-wgt'"))

        l_create = ttk.Checkbutton(self.lfr_ctpl, text="     create tpl", variable=self.cb_createtpl, command=self.Update_ctpl)
        self.cb_createtpl.set(self.configs.get('createtpl', 0))
        l_create.grid(row=0, column=0, sticky="nw", padx=(xy0,x1), pady=y1, columnspan=2)
        Help_Box(widget = l_create, text = text_from_file("'-createtpl'"))

        l_dat = ttk.Checkbutton(lfr_out, text="  .dat", variable=self.cb_format[0])
        l_dat.grid(row=1, column=1, sticky="nw", padx=(xy0,5), pady=y1)
        self.cb_format[0].set(1)

        l_fits = ttk.Checkbutton(lfr_out, text="  .fits (astropy)", variable=self.cb_format[1])
        l_fits.grid(row=1, column=2, sticky="nw", padx=5, pady=y1)

        l_cpl = ttk.Checkbutton(lfr_out, text="  .fits (CPL)", variable=self.cb_format[2])
        l_cpl.grid(row=1, column=3, sticky="nw", padx=5, pady=y1)

        self.combo_ip = ttk.Combobox(lfr_model, values=[*IPs])	# IPs from model.py
        self.combo_ip.set(self.configs.get('ip','g'))
        self.combo_ip.grid(row=0, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

        self.combo_tell = ttk.Combobox(self.lfr_tell, values=['', 'add', 'add2', 'mask', 'sig'], width=14)
        self.combo_tell.set(self.configs.get('telluric','mask'))
        self.combo_tell.grid(row=0, column=1, sticky="nw", padx=(x1), pady=y1, columnspan=2)
        self.combo_tell.bind('<<ComboboxSelected>>', lambda event: self.Update_tell())

    def fr_plot(self):
        # Frame for parameter selection plotting
        fr3 = Frame(self.master, bg=bg_frame, bd=2, relief='groove')
        fr3.grid(row=1, column = 1, sticky="news", padx=(0, xy0), pady=0)
        fr3.grid_propagate(False)
        fr3.grid_columnconfigure(0, weight=1)

        lfr_plot1 = LabelFrame(fr3, text="Data Analysis", bg=bg_frame, bd=2)
        lfr_plot1.grid(row=1, column=0, sticky="news", padx=(10,10), pady=y1, ipady=5)

        lfr_plot2 = LabelFrame(fr3, text="Plot fitted chunks", bg=bg_frame, bd=2)
        lfr_plot2.grid(row=2, column=0, sticky="news", padx=(10,10), pady=y1, ipady=5)

        lfr_plot1.grid_columnconfigure(0, weight=1)
        lfr_plot1.grid_columnconfigure(1, weight=1)
        lfr_plot1.grid_columnconfigure(2, weight=1)

        lfr_plot2.grid_columnconfigure(0, weight=1)
        lfr_plot2.grid_columnconfigure(1, weight=1)

        # Label
        l_plot = ttk.Label(fr3, text='Options plotting data', font=(font_type, font_size, 'bold'))
        l_plot.grid(row=0, column=0, sticky="nw", padx=(xy0,0), pady=(10,y1))

        # Entry
        self.e_lookpar = Entry(lfr_plot1)
        self.e_lookpar.grid(row=5, column=1, sticky="news", padx=(10,xy0), pady=y1, columnspan=2)

        self.e_lookfast = Entry(lfr_plot2)
        self.e_lookfast.grid(row=0, column=1, sticky="news", padx=(0,xy0), pady=y1, columnspan=2)

        self.e_look = Entry(lfr_plot2)
        self.e_look.grid(row=1, column=1, sticky="news", padx=(0,xy0), pady=y1, columnspan=2)

        # Checkboxes
        self.cb_demo = [IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar()]
        self.Checkbox_list(lfr_plot1, ["raw data", "plot IP", "stellar tpl", "forward model", "fit continuum", "wavelength solution", "fit vguess"], self.cb_demo)

        self.cb_lookguess = IntVar()
        self.cb_lookres = IntVar()
        self.cb_infoprec = IntVar()
        self.cb_lookfast = IntVar()
        self.cb_lookpar = IntVar()
        self.cb_look = IntVar()

        l_lguess = ttk.Checkbutton(lfr_plot1, text="     lookguess", variable=self.cb_lookguess)
        l_lguess.grid(row=3, column=2, sticky="nw", padx=(xy0,x1), pady=y1)
        Help_Box(widget = l_lguess, text = "Plot model using start guess of parameters.")

        l_res = ttk.Checkbutton(lfr_plot1, text="     lookres", variable=self.cb_lookres)
        l_res.grid(row=4, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
        Help_Box(widget = l_res, text = text_from_file("'-lookres'"))

        l_infop = ttk.Checkbutton(lfr_plot1, text="     infoprec", variable=self.cb_infoprec)
        l_infop.grid(row=4, column=2, sticky="nw", padx=(xy0,x1), pady=y1)
        Help_Box(widget = l_infop, text = text_from_file("'-infoprec'"))

        l_lpar = ttk.Checkbutton(lfr_plot1, text="     lookpar:", variable=self.cb_lookpar)
        l_lpar.grid(row=5, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
        Help_Box(widget = l_lpar, text = text_from_file("'-lookpar'"))

        l_lfast = ttk.Checkbutton(lfr_plot2, text="     lookfast:", variable=self.cb_lookfast)
        l_lfast.grid(row=0, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
        Help_Box(widget = l_lfast, text = text_from_file("'-lookfast'"))
        self.cb_lookfast.set(1)

        l_look = ttk.Checkbutton(lfr_plot2, text="     look:", variable=self.cb_look)
        l_look.grid(row=1, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
        Help_Box(widget = l_look, text = text_from_file("'-look'"))

    def bt_file(self, e_file):
        file = askopenfilename()
        if file:
            e_file.delete(0, END)
            e_file.insert(0, file)

    def Update_inst(self):
        Inst = importlib.import_module('inst.inst_'+self.combo_inst.get())
        FTS = Inst.FTS
     #   default=viperdir + FTS.__defaults__[0]
        default=FTS.__defaults__[0]
        self.e_cell.delete(0, END)
        self.e_cell.insert(0, default)

        self.e_iset.delete(0, END)
        self.e_iset.insert(0, getattr(Inst, 'iset', ':'))

        # molecule selection for tellurics
        #self.molec = list(Inst.atmall.keys())

        if str(self.combo_tell.get()) in ('add', 'add2'):
            self.Telluric()

    def Update_ctpl(self):
        if hasattr(self, 'l_kapctpl'): self.l_kapctpl.destroy()
        if hasattr(self, 'e_kapctpl'): self.e_kapctpl.destroy()
        if hasattr(self, 'l_lookctpl'): self.l_lookctpl.destroy()

        if self.cb_createtpl.get():
            self.l_kapctpl = ttk.Label(self.lfr_ctpl,text='kapsig_ctpl:')
            self.l_kapctpl.grid(row=2, column=0, sticky="nw", padx=(xy0,0), pady=y1)
            Help_Box(widget = self.l_kapctpl, text = text_from_file("'-kapsig_ctpl'"))

            self.e_kapctpl = Entry(self.lfr_ctpl, width=9)
            self.e_kapctpl.insert(0, self.configs.get('kapsig_cptl','0.6'))
            self.e_kapctpl.grid(row=2, column=1, sticky="nw", padx=(x1,xy0), pady=y1)#, columnspan=2)

            self.l_lookctpl = ttk.Checkbutton(self.lfr_ctpl, text="     lookctpl", variable=self.cb_lookctpl)
            self.l_lookctpl.grid(row=1, column=0, sticky="nw", padx=(xy0,x1), pady=y1)#, columnspan=2)
            Help_Box(widget = self.l_lookctpl, text = text_from_file("'-lookctpl'"))
            self.cb_lookctpl.set(1)

    def Update_tell(self):
        if hasattr(self, 'l_tsig'): self.l_tsig.destroy()
        if hasattr(self, 'e_tsig'): self.e_tsig.destroy()
        if hasattr(self, 'l_tshift'): self.l_tshift.destroy()

        if str(self.combo_tell.get()) in ('add', 'add2', 'sig'):
            self.l_tsig = ttk.Label(self.lfr_tell,text='tsig:')
            self.l_tsig.grid(row=1, column=0, sticky="nw", padx=(xy0,0), pady=y1)
            Help_Box(widget = self.l_tsig, text = text_from_file("'-tsig'"))

            self.e_tsig = Entry(self.lfr_tell, width=8)
            self.e_tsig.insert(0, self.configs.get('tsig','1'))
            self.e_tsig.grid(row=1, column=1, sticky="nw", padx=(x1,xy0), pady=y1, columnspan=2)

        if str(self.combo_tell.get()) in ('add', 'add2'):
            self.l_tshift = ttk.Checkbutton(self.lfr_tell, text="     tell shift", variable=self.cb_tellshift)
            self.cb_tellshift.set(self.configs.get('tellshift', 0))
            self.l_tshift.grid(row=2, column=0, sticky="nw", padx=(xy0,x1), pady=y1, columnspan=3)
            Help_Box(widget = self.l_tshift, text = text_from_file("'-tellshift'"))

            self.Telluric(show=1)
        else:
            self.Telluric(show=0)

    def Telluric(self, show=1):
        if self.cb_atm != []:
            # delete old checkboxes, label
            for i in self.cb_atm: i.destroy()
            self.cb_atm.clear()
            self.cbv_atm.clear()
            self.l_molec.destroy()

        if show:
            if str(self.combo_inst.get()) in ('TLS', 'CES', 'OES', 'KECK', 'UVES', 'McDonald'):
                self.l_molec = ttk.Label(self.lfr_tell, text='Optical molecules:')
                self.l_molec.grid(row=3, column=0, sticky="nw", padx=(xy0,0), pady=y1, columnspan=6)
                self.molec = ['H2O', 'O2']
            elif str(self.combo_inst.get()) in ('CRIRES', 'cplCRIRES'):
                self.l_molec = ttk.Label(self.lfr_tell, text='Near Infrared molecules:')
                self.l_molec.grid(row=3, column=0, sticky="nw", padx=(xy0,0), pady=y1, columnspan=6)
                self.molec = ['H2O', 'CH4', 'N2O', 'CO2', 'CO']

            for m, mol in enumerate(self.molec):
                yi, xi =  divmod(m, 3)
                self.cbv_atm.append(IntVar())
                c = ttk.Checkbutton(self.lfr_tell, text=" "+mol, variable=self.cbv_atm[m])
                c.grid(row=4+yi, column=xi, sticky="news", padx=(20,0), pady=(y1,0))
                self.cb_atm.append(c)
                self.cbv_atm[m].set(1)


    def Label_list(self, frame, lab_list):
        for l, labt in enumerate (lab_list):
           lab = ttk.Label(frame, text=str(labt)+':')
           lab.grid(row=l, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
           Help_Box(widget = lab, text = text_from_file("'-"+str(labt)+"'"))

    def Checkbox_list(self, frame, cb_list, cbv):
        for c, cbt in enumerate (cb_list):
            co, ro =  divmod(c, 4)
            if co == 1: co = 2
            l_raw = ttk.Checkbutton(frame, text="     "+str(cbt), variable=cbv[c])
            l_raw.grid(row=ro, column=co, sticky="nw", padx=(xy0,x1), pady=y1)

    def bt_start(self):
        print('---Starting viper.py with choosen parameters---')

        ar_demo = np.zeros(len(self.cb_demo))
        for d, dem in enumerate(self.cb_demo):
            ar_demo[d] = dem.get()
        x = 2 * np.ones(len(self.cb_demo))
        demo = np.poly1d(ar_demo[::-1])(x)[0]

        str_arg = self.e_dat.get()+"' " + self.e_tpl.get()
        str_arg += " -inst " + self.combo_inst.get()
        str_arg += " -nset " + str(self.e_nset.get()) 
        str_arg += " -oset "+self.e_oset.get()
        str_arg += " -chunks " + self.e_ch.get()
        str_arg += " -vcut " + self.e_vcut.get() 
        str_arg += " -iset "+self.e_iset.get()
        str_arg += " -ip " + self.combo_ip.get() 
        str_arg += " -iphs " + self.e_iphs.get()
        str_arg += " -deg_norm " + self.e_deg_norm.get() 
        str_arg += " -deg_wave " + self.e_deg_wave.get() 
        str_arg += " -deg_bkg " + self.e_deg_bkg.get()
        str_arg += " -rv_guess " + self.e_vg.get()
        str_arg += " -demo " + str(int(demo))

        if self.e_targ.get(): str_arg += " -targ " + str(self.e_targ.get())
        if self.e_tag.get(): str_arg += " -tag " + str(self.e_tag.get())
        if self.e_kapsig.get(): str_arg += " -kapsig " + str(self.e_kapsig.get()) 
        if self.cb_wgt.get(): str_arg += " -wgt "
        if self.e_overs.get(): str_arg += " -oversampling " + str(self.e_overs.get())
        if self.cb_lookpar.get(): str_arg += " -lookpar " + self.e_lookpar.get()
        if self.cb_lookguess.get(): str_arg += " -lookguess "
        if self.cb_lookres.get(): str_arg += " -lookres "
        if self.cb_infoprec.get(): str_arg += " -infoprec "
        if self.cb_look.get(): str_arg += " -look " + self.e_look.get()
        if self.cb_lookfast.get(): str_arg += " -lookfast " + self.e_lookfast.get()

        if self.e_flag.get() != '' and self.cb_flagfile.get():
             str_arg += " -flagfile "+ str(self.e_flag.get())

        if self.cb_cell.get(): 
            # str_arg += " -fts " + viperdir + self.e_cell.get()
             str_arg += " -fts " + self.e_cell.get()
        else:
             str_arg += " -nocell " 

        if self.combo_tell.get():
            str_arg += " -telluric " + str(self.combo_tell.get())
        if str(self.combo_tell.get()) in ('add', 'add2', 'sig'):
            str_arg += " -tsig " + self.e_tsig.get() 
        if str(self.combo_tell.get()) in ('add', 'add2'):
            if self.cb_tellshift.get(): str_arg += " -tellshift "
            str_arg += " -molec "
            for m, atm in enumerate(self.cbv_atm):        
                if atm.get():
                    str_arg += str(self.molec[m])+" "

        if self.cb_createtpl.get(): 
             str_arg += " -createtpl "
             if self.cb_lookctpl.get(): str_arg += " -lookctpl "
             str_arg += " -kapsig_ctpl " + str(self.e_kapctpl.get()) 

        oformat = np.array([f.get() for f in self.cb_format])
        str_arg += " -output_format " + " ".join((np.array(["dat", "fits", "cpl"])[oformat==1]))

        #self.cb_format

        self.e_run.delete('0.0', END)
        self.e_run.insert(INSERT,"python3 viper.py '"+str_arg)
        self.e_run.update()

        print(str_arg)

        os.system("python3 "+viperdir+"viper.py '" + str_arg)

        print('---Finished viper.py---')


def main():
    win = Tk()
    win.title('Gui VIPER')
    win.geometry("{}x{}".format(win_width, win_high))
    win.configure(background=bg_color)

    win.grid_propagate(False)
    win.columnconfigure(0, weight=1, uniform=1)
    win.columnconfigure(1, weight=1, uniform=1)
    win.rowconfigure(0, weight=5, uniform=1, minsize=188)
    win.rowconfigure(1, weight=5, uniform=1, minsize=338)
    win.rowconfigure(2, weight=1, uniform=1)
    win.rowconfigure(3, weight=5, uniform=1) #, minsize=130)
    win.rowconfigure(4, weight=2, uniform=1)

    s = ttk.Style()
    s.configure("TCheckbutton", background=bg_frame, bd=0, highlightthickness=0)  
    s.configure('TRadiobutton', background=bg_frame)   
    s.configure("bold.TButton", padding=2, background=bg_button, font=(font_type,font_size,'bold'), borderwidth=2, width=15)
    s.configure("nbold.TButton", padding=2, background=bg_button, font=(font_type,font_size-1,''), borderwidth=2, width=15)     
    s.configure('TFrame', background=bg_frame)
    s.configure('TLabel', background=bg_frame)

    configs = {}    
    # Using config.ini file for reading in parameter values
    if len(sys.argv) == 3 and sys.argv[1].endswith('.ini'):
            
        config_file = sys.argv[1] #'config_viper.ini'
        config_sect = sys.argv[2] #'CRIRES'

        config = configparser.ConfigParser()
        config.read(config_file)
        if str(config_sect) not in config.sections():
            print('Section not found in parser file. Using default values')
            config_sect = 'DEFAULT'
            
        configs = dict(config[str(config_sect)])  
        if 'kapsig' in configs.keys():
            configs['kapsig'] = [float(i) for i in (configs['kapsig']).split(' ')]      

    GUI_viper(win, configs)
    
    win.mainloop()

if __name__ == '__main__':
    main()
