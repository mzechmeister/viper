#! /usr/bin/env python3
# -*- coding: iso-8859-1 -*-

# GUI to start viper.py, showing most important options

#import gplot, gnuplot.funcutils
import os
import subprocess
from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter import ttk
from gplot import *
import numpy as np
import importlib
from tkinter.scrolledtext import ScrolledText
from hbox import Help_Box

viperdir = os.path.dirname(os.path.realpath(__file__)) + os.sep

gplot.colors('classic')
gplot2 = Gplot()
g = Gplot()


###### FUNCIONS BUTTONS ######

def bt_start():
    print('---Starting viper.py with choosen parameters---')

    ar_demo = np.zeros(len(cb_demo))
    for d, dem in enumerate(cb_demo):
        ar_demo[d] = dem.get()
    x = 2 * np.ones(len(cb_demo))
    demo = np.poly1d(ar_demo[::-1])(x)[0]

    str_arg = e_dat.get()+"' "+e_tpl.get()+" -inst "+combo_inst.get()+" -deg_norm "+e_deg_norm.get()+" -deg_wave "+e_deg_wave.get()+" -deg_bkg "+e_deg_bkg.get()+" -nset "+str(e_nset.get())+" -oset "+e_oset.get()+" -chunks "+e_ch.get()+" -demo "+str(int(demo))+" -rv_guess "+e_vg.get() +" -ip "+combo_ip.get() +" -iphs "+e_iphs.get() +" -tsig "+ e_tsig.get() +" -vcut "+e_vcut.get() 

 #   if combo_stepRV.get():
  #      str_arg += " -stepRV " + str(combo_stepRV.get())
    if e_kapsig.get():
        str_arg += " -kapsig "+ str(e_kapsig.get()) 
    if cb_lookpar.get():
        str_arg += " -lookpar " + e_lookpar.get()
    if cb_look.get():
        str_arg += " -look " + e_look.get()
    if cb_lookfast.get():
        str_arg += " -lookfast " + e_lookfast.get()
    if cb_lookguess.get():
        str_arg += " -lookguess "
    if cb_lookres.get():
        str_arg += " -lookres "
    if cb_infoprec.get():
        str_arg += " -infoprec "
    if e_overs.get():
        str_arg += " -oversampling "+ str(e_overs.get())
    if e_targ.get():
        str_arg += " -targ " + str(e_targ.get())
    if e_tag.get():
        str_arg += " -tag " + str(e_tag.get())
    if e_flag.get() != '' and cb_flagfile.get():
         str_arg += " -flagfile "+ str(e_flag.get())
  #  if cb_nocell.get():
   #     str_arg += " -nocell "
    if cb_createtpl.get():
        str_arg += " -createtpl "
    if cb_tellshift.get():
        str_arg += " -tellshift "
    if cb_wgt.get():
        str_arg += " -wgt "
    if cb_cell.get():
        str_arg += " -fts " + viperdir + e_cell.get()
    else:
        str_arg += " -nocell " 
    if combo_tell.get():
        str_arg += " -telluric " + str(combo_tell.get())
        str_arg += " -molec "
        if combo_inst.get() in ('TLS', 'CES', 'OES', 'KECK'):
            cb_atm2, molec2 = cb_atm[:2], molec[:2]
        elif combo_inst.get() in ('CRIRES', 'cplCRIRES'):
            cb_atm2, molec2 = cb_atm[2:], molec[2:]
        for m, atm in enumerate(cb_atm2):        
            if atm.get():
                str_arg += str(molec2[m])+" "

    e_run.delete('0.0', END)
    e_run.insert(INSERT,"python3 viper.py "+str_arg)
    e_run.update()

    print(str_arg)

    os.system("python3 "+viperdir+"viper.py '" + str_arg)

   # doesn't work as subprocess
   # os.system("python3 viper.py 'data/CRIRES/Ross619/210*' data/CRIRES/model/Ross619_19.fits -inst CRIRES -deg_norm 2 -deg_wave 2 -nset :2 -oset 10")
  #  subprocess.Popen(["./viper.py","'data/CRIRES/Ross619/210*'", "data/CRIRES/model/Ross619_19.fits -inst CRIRES -deg_norm 2 -deg_wave 2 -nset :2 -oset 10 -test 1"])

    print('---Finished viper.py---')

def text_from_file(text):
     with open(viperdir+"viper.py") as search:
         for line in search:
             line = line.rstrip()  
             if text in line:
                  infotext = (line.split('help')[1]).split("'")[1]
     return infotext

def bt_file(e_file):
    file = askopenfilename(master=win)
    if file:
        e_file.delete(0, END)
        e_file.insert(0, file)

def bt_exit():
    exit()

def new(event):
    win_width = win.winfo_width()
    win_height = event.height
    fr1.config(width=win_width-40)
    fr2.config(width=(win_width-46)/2)
    fr3.config(width=(win_width-46)/2)

def FTS_default(*args):
    Inst = importlib.import_module('inst.inst_'+combo_inst.get())
    FTS = Inst.FTS
 #   default=viperdir + FTS.__defaults__[0]
    default=FTS.__defaults__[0]
    e_cell.delete(0, END)
    e_cell.insert(0, default)

###### SETTING ######

font_type = 'Arial'
font_size = 12
bg_frame = '#f1f1f1'    # bg color big frame
bg_color = '#e6e1e1'    # bg color small frames
bg_button = "#fdfdfd"
#bg_frame = 'SlateGray1'	# bg color small frames
#bg_color = 'SkyBlue3'		# bg color big frame
#bg_button = 'alice blue'	# azure

win_width = 1060     # width of GUI window	
win_high = 780      # height of GUI window	
xy0 = 20
y1 = 3
x1 = 10

win = Tk()
win.title('Gui VIPER')
win.geometry("{}x{}".format(win_width, win_high))
win.configure(background=bg_color)


###### FRAMES #######
win.bind("<Configure>", new)
win.columnconfigure(0, weight=1)
win.columnconfigure(2, weight=1)
win.rowconfigure(3, weight=1)

# options: solid, sunken, ridge, groove, raised
# groove, ridge need bd=2 or more

# Frame for data input
fr1 = Frame(win, height=168, width=win_width-40, bg=bg_frame, bd=2, relief='groove')
fr1.grid(row=0, column = 1, sticky="nw",padx=20,pady=(xy0,6), columnspan=2)
fr1.grid_propagate(False)
fr1.grid_columnconfigure(2, weight=1)
fr1.grid_columnconfigure(4, weight=1)
fr1.grid_rowconfigure(0, weight=1)
fr1.grid_rowconfigure(5, weight=1)

# Frame for parameter selection data reduction
fr2 = Frame(master=win, height=568, width=(win_width-46)/2, bg=bg_frame, bd=2, relief='groove')
fr2.grid(row=1, column = 1, sticky="nw",padx=(xy0,6),pady=0, rowspan=4)
fr2.grid_propagate(False)
fr2.grid_columnconfigure(0, weight=1)
fr2.grid_columnconfigure(1, weight=1)

# Frame for parameter selection plotting
fr3 = Frame(master=win, height=338, width=(win_width-46)/2, bg=bg_frame, bd=2, relief='groove')
fr3.grid(row=1, column = 2, sticky="nw",padx=(0, xy0),pady=0)
fr3.grid_propagate(False)
fr3.grid_columnconfigure(0, weight=1)

# Sub-Frames of fr2
lfr_data = LabelFrame(fr2, text="Data", bg=bg_frame, bd=2)#, width=win_width-40)
lfr_data.grid(row=1, column=0, sticky="news", padx=(10,0), pady=y1, ipady=5)

lfr_model = LabelFrame(fr2, text="Model", bg=bg_frame, bd=2)#, width=win_width-40)
lfr_model.grid(row=1, column=1, sticky="news", padx=(10,10), pady=y1, ipady=5)

lfr_tpl = LabelFrame(fr2, text="Template", bg=bg_frame, bd=2)#, width=win_width-40)
lfr_tpl.grid(row=2, column=0, sticky="news", padx=(10,0), pady=y1, ipady=5)

lfr_stat = LabelFrame(fr2, text="Fit Parameters", bg=bg_frame, bd=2)#, width=win_width-40)
lfr_stat.grid(row=2, column=1, sticky="news", padx=(10,10), pady=y1, ipady=5)

lfr_tell = LabelFrame(fr2, text="Tellurics", bg=bg_frame, bd=2)#, width=win_width-40)
lfr_tell.grid(row=3, column=0, sticky="news", padx=(10,10), pady=y1, ipady=5, columnspan=2)

for lfr in [lfr_data, lfr_tpl, lfr_model, lfr_stat, lfr_tell]:
    lfr.grid_columnconfigure(1, weight=1)

for c in range(0,5,1):
    lfr_tell.grid_columnconfigure(c, weight=1)

# Sub-Frames of fr3
lfr_plot1 = LabelFrame(fr3, text="Data Analysis", bg=bg_frame, bd=2)#, width=win_width-40)
lfr_plot1.grid(row=1, column=0, sticky="news", padx=(10,10), pady=y1, ipady=5)

lfr_plot2 = LabelFrame(fr3, text="Plot fitted chunks", bg=bg_frame, bd=2)#, width=win_width-40)
lfr_plot2.grid(row=2, column=0, sticky="news", padx=(10,10), pady=y1, ipady=5)

lfr_plot1.grid_columnconfigure(0, weight=1)
lfr_plot1.grid_columnconfigure(1, weight=1)
lfr_plot1.grid_columnconfigure(2, weight=1)

lfr_plot2.grid_columnconfigure(0, weight=1)
lfr_plot2.grid_columnconfigure(1, weight=1)

###### BUTTONS ######

ttk.Style().configure("TButton", padding=2,   background=bg_button, font=(font_type,font_size,'bold'), borderwidth =2)

b_exit = ttk.Button(master=win, text='EXIT', command = bt_exit)
b_exit.grid(row=4, column=2, sticky="se", padx=(0,xy0), pady = (10,20))

b_go = ttk.Button(master=win, text='Start', command = bt_start)
b_go.grid(row=4, column=2, sticky="se", padx=(0,130), pady = (10,20))

b_dat = Button(fr1, text='Search data file', command = lambda: bt_file(e_dat), background=bg_button, width=15)
b_dat.grid(row=1, column=7, sticky="nw", padx=xy0)

b_tpl = Button(fr1,text='Search tpl file', command = lambda: bt_file(e_tpl), background=bg_button, width=15)
b_tpl.grid(row=2, column=7, sticky="nw", padx=xy0)

b_cell = Button(fr1,text='Search Cell file', command = lambda: bt_file(e_cell), background=bg_button, width=15)
b_cell.grid(row=3, column=7, sticky="nw", padx=xy0)

b_flag = Button(fr1,text='Search flag file', command = lambda: bt_file(e_flag), background=bg_button, width=15)
b_flag.grid(row=4, column=7, sticky="nw", padx=xy0)

###### ENTRIES ######

e_dat = Entry(fr1, width = 200)
e_dat.insert(0, 'data/TLS/hd189733/*')
e_dat.grid(row=1, column=1, sticky="nw", padx=x1, pady=y1, columnspan=6)

e_tpl = Entry(fr1, width = 200)
e_tpl.insert(0, 'data/TLS/Deconv/HARPS*fits')
e_tpl.grid(row=2, column=1, sticky="nw", padx=x1, pady=y1, columnspan=6)

e_cell = Entry(fr1, width = 200)
e_cell.grid(row=3, column=1, sticky="nw", padx=x1, pady=y1, columnspan=6)

e_flag = Entry(fr1, width = 200)
e_flag.insert(0, '')
e_flag.grid(row=4, column=1, sticky="nw", padx=x1, pady=y1, columnspan=6)

e_targ = Entry(fr1, width=70)
e_targ.grid(row=5, column=4, sticky="nw", padx=x1, pady=y1)

e_tag = Entry(fr1, width=35)
e_tag.grid(row=5, column=6, sticky="nw", padx=xy0, pady=y1, columnspan=2)

e_nset = Entry(lfr_data)#, width=15)
e_nset.insert(0, ':1')
e_nset.grid(row=0, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

e_oset = Entry(lfr_data)
e_oset.insert(0, '20')
e_oset.grid(row=1, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

e_ch = Entry(lfr_data)
e_ch.insert(0, '1')
e_ch.grid(row=2, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

e_vcut = Entry(lfr_data)
e_vcut.insert(0, '100')
e_vcut.grid(row=3, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

e_iphs = Entry(lfr_model)
e_iphs.insert(0, '50')
e_iphs.grid(row=1, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

e_deg_norm = Entry(lfr_model)
e_deg_norm.insert(0, '3')
e_deg_norm.grid(row=2, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

e_deg_wave = Entry(lfr_model)
e_deg_wave.insert(0, '3')
e_deg_wave.grid(row=3, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

e_deg_bkg = Entry(lfr_model)
e_deg_bkg.insert(0, '1')
e_deg_bkg.grid(row=4, column=1, sticky="nw", padx=(x1,xy0), pady=(y1,0))

e_vg = Entry(lfr_tpl)
e_vg.insert(0, '1')
e_vg.grid(row=0, column=1, sticky="nw", padx=(x1,xy0), pady=(0,0))

e_overs = Entry(lfr_tpl)
e_overs.insert(0, '1')
e_overs.grid(row=1, column=1, sticky="nw", padx=(x1,xy0), pady=0)

e_kapsig = Entry(lfr_stat)
e_kapsig.insert(0, '4.5')
e_kapsig.grid(row=0, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

e_tsig = Entry(lfr_tell, width=8)
e_tsig.insert(0, '1')
e_tsig.grid(row=1, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

e_lookpar = Entry(lfr_plot1)
e_lookpar.grid(row=6, column=1, sticky="news", padx=(10,xy0), pady=y1, columnspan=2)

e_lookfast = Entry(lfr_plot2)
e_lookfast.grid(row=8, column=1, sticky="news", padx=(0,xy0), pady=y1, columnspan=2)

e_look = Entry(lfr_plot2)
e_look.grid(row=9, column=1, sticky="news", padx=(0,xy0), pady=y1, columnspan=2)

e_run = ScrolledText(win, background=bg_frame)
e_run.grid(row=3, column = 2, sticky="news",padx=(0, xy0),pady=(0,5))

###### COMBOBOXES ######

combo_inst = ttk.Combobox(fr1, values=['TLS', 'CRIRES','cplCRIRES', 'CES', 'KECK', 'UVES', 'OES'], width=10)
combo_inst.set('TLS')
combo_inst.grid(row=5, column=1, sticky="nw", padx=x1, pady=y1)
combo_inst.bind('<<ComboboxSelected>>', lambda event: FTS_default())

combo_ip = ttk.Combobox(lfr_model, values=['g', 'ag', 'agr', 'bg', 'mcg', 'mg', 'sg', 'bnd'])
combo_ip.set('g')
combo_ip.grid(row=0, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

#combo_stepRV = ttk.Combobox(fr2, values=['', 'a', 'm'])
#combo_stepRV.set('')
#combo_stepRV.grid(row=6, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

combo_tell = ttk.Combobox(lfr_tell, values=['', 'add', 'add2', 'mask', 'sig'], width=6)
combo_tell.set('')
combo_tell.grid(row=0, column=1, sticky="nw", padx=(x1), pady=y1, columnspan=2)

###### CHECKBOXES ######

ttk.Style().configure("TCheckbutton", background=bg_frame, bd=0, highlightthickness=0)

cb_cell = IntVar()
l_cell = ttk.Checkbutton(fr1, text="   Cell file:", variable=cb_cell)
l_cell.grid(row=3, column=0, sticky="nw", padx=x1, pady=y1)
Help_Box(widget = l_cell, text = text_from_file("'-fts'") + " If checkbox is unset, no use of cell FTS for the modelling.")
cb_cell.set(1) 

cb_flagfile = IntVar()
l_flag = ttk.Checkbutton(fr1, text="   flag file:", variable=cb_flagfile)
l_flag.grid(row=4, column=0, sticky="nw", padx=x1, pady=y1)
Help_Box(widget = l_flag, text = text_from_file("'-flagfile'"))

cb_wgt = IntVar()
l_wei = ttk.Checkbutton(lfr_stat, text="     weighted error", variable=cb_wgt)
l_wei.grid(row=1, column=0, sticky="nw", padx=(xy0,x1), pady=y1, columnspan=2)
Help_Box(widget = l_wei, text = text_from_file("'-wgt'"))

#cb_nocell = IntVar()
#ttk.Checkbutton(fr2, text="     no cell", variable=cb_nocell).grid(row=x1, column=2, sticky="nw", padx=(xy0,x1), pady=y1)
cb_createtpl = IntVar()
l_create = ttk.Checkbutton(lfr_tpl, text="     create tpl", variable=cb_createtpl)
l_create.grid(row=2, column=0, sticky="nw", padx=(xy0,x1), pady=y1, columnspan=2)
cb_tellshift = IntVar()
Help_Box(widget = l_create, text = text_from_file("'-createtpl'"))

l_tshift = ttk.Checkbutton(lfr_tell, text="     tell shift", variable=cb_tellshift)
l_tshift.grid(row=2, column=0, sticky="nw", padx=(xy0,x1), pady=y1, columnspan=3)
Help_Box(widget = l_tshift, text = text_from_file("'-tellshift'"))

cb_atm = [IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar()]
ttk.Checkbutton(lfr_tell, text="   H2O", variable=cb_atm[0]).grid(row=4, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(lfr_tell, text="   O2", variable=cb_atm[1]).grid(row=4, column=1, sticky="nw", padx=(xy0,x1), pady=y1)

ttk.Checkbutton(lfr_tell, text="   H2O", variable=cb_atm[2]).grid(row=6, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(lfr_tell, text="   CH4", variable=cb_atm[3]).grid(row=6, column=1, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(lfr_tell, text="   CO", variable=cb_atm[4]).grid(row=6, column=2, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(lfr_tell, text="   CO2", variable=cb_atm[5]).grid(row=6, column=3, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(lfr_tell, text="   N2O", variable=cb_atm[6]).grid(row=6, column=4, sticky="nw", padx=(xy0,x1), pady=y1)

molec = ['H2O', 'O2', 'H2O', 'CH4', 'CO', 'CO2', 'N2O']
for mol in cb_atm:
    mol.set(1)

# plotting options
cb_demo = [IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar()]

l_raw = ttk.Checkbutton(lfr_plot1, text="     raw data", variable=cb_demo[0])
l_raw.grid(row=1, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_raw, text = "Plot the raw input data, template and cell.")
l_pip = ttk.Checkbutton(lfr_plot1, text="     plot IP", variable=cb_demo[1])
l_pip.grid(row=2, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_pip, text = "Plot start guess for IP.")
l_stpl = ttk.Checkbutton(lfr_plot1, text="     stellar tpl", variable=cb_demo[2])
l_stpl.grid(row=3, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_stpl, text = "Plot model for stellar template.")
ttk.Checkbutton(lfr_plot1, text="     forward model", variable=cb_demo[3]).grid(row=4, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(lfr_plot1, text="     fit continuum", variable=cb_demo[4]).grid(row=1, column=2, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(lfr_plot1, text="     wavelength solution", variable=cb_demo[5]).grid(row=2, column=2, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(lfr_plot1, text="     fit vguess", variable=cb_demo[6]).grid(row=3, column=2, sticky="nw", padx=(xy0,x1), pady=y1)
cb_lookguess = IntVar()
l_lguess = ttk.Checkbutton(lfr_plot1, text="     lookguess", variable=cb_lookguess)
l_lguess.grid(row=4, column=2, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_lguess, text = "Plot model using start guess of parameters.")
cb_lookres = IntVar()
l_res = ttk.Checkbutton(lfr_plot1, text="     lookres", variable=cb_lookres)
l_res.grid(row=5, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_res, text = text_from_file("'-lookres'"))
cb_infoprec = IntVar()
l_infop = ttk.Checkbutton(lfr_plot1, text="     infoprec", variable=cb_infoprec)
l_infop.grid(row=5, column=2, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_infop, text = text_from_file("'-infoprec'"))
cb_lookfast = IntVar()
cb_lookpar = IntVar()
l_lpar = ttk.Checkbutton(lfr_plot1, text="     lookpar:", variable=cb_lookpar)
l_lpar.grid(row=6, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_lpar, text = text_from_file("'-lookpar'"))

l_lfast = ttk.Checkbutton(lfr_plot2, text="     lookfast:", variable=cb_lookfast)
l_lfast.grid(row=8, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_lfast, text = text_from_file("'-lookfast'"))
cb_lookfast.set(1)
cb_look = IntVar()
l_look = ttk.Checkbutton(lfr_plot2, text="     look:", variable=cb_look)
l_look.grid(row=9, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_look, text = text_from_file("'-look'"))

###### LABELS ######

Label(fr1, text='data files', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=1, column=0, sticky="nw", padx=x1, pady=y1)
Label(fr1, text='template file', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=2, column=0, sticky="nw", padx=x1, pady=y1)

l_inst = Label(fr1, text='spectrograph:', background=bg_frame)
l_inst.grid(row=5, column=0, sticky="nw", padx=x1, pady=y1)
Help_Box(widget = l_inst, text = text_from_file("'-inst'"))

l_targ = Label(fr1, text='targ:', background=bg_frame)
l_targ.grid(row=5, column=3, sticky="nw", padx=xy0, pady=y1)
Help_Box(widget = l_targ, text = text_from_file("'-targ'"))
l_tag = Label(fr1, text='tag:', background=bg_frame)
l_tag.grid(row=5, column=5, sticky="nw", padx=xy0, pady=y1)
Help_Box(widget = l_tag, text = text_from_file("'-tag'"))

Label(fr2, text='Options data reduction', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=0, column=0, sticky="nw", padx=(xy0,0), pady=(10,y1), columnspan=3)
Label(fr3, text='Options plotting data', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=0, column=0, sticky="nw", padx=(xy0,0), pady=(10,y1))

l_nset = Label(lfr_data, text='nset:', background=bg_frame)
l_nset.grid(row=0, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_nset, text = text_from_file("'-nset'"))
l_oset = Label(lfr_data, text='oset:', background=bg_frame)
l_oset.grid(row=1, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_oset, text = text_from_file("'-oset'"))
l_ch = Label(lfr_data, text='chunks:', background=bg_frame)
l_ch.grid(row=2, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_ch, text = text_from_file("'-chunks'"))
l_vcut = Label(lfr_data,text='vcut:', background=bg_frame)
l_vcut.grid(row=3, column=0, sticky="nw", padx=(xy0,x1), pady=(y1))
Help_Box(widget = l_vcut, text = text_from_file("'-vcut'"))

l_ip = Label(lfr_model,text='IP:', background=bg_frame)
l_ip.grid(row=0, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_ip, text = text_from_file("'-ip'"))
l_iphs = Label(lfr_model,text='iphs:', background=bg_frame)
l_iphs.grid(row=1, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_iphs, text = text_from_file("'-iphs'"))
l_norm = Label(lfr_model,text='deg_norm:', background=bg_frame)
l_norm.grid(row=2, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_norm, text = text_from_file("'-deg_norm'"))
l_wave = Label(lfr_model,text='deg_wave:', background=bg_frame)
l_wave.grid(row=3, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_wave, text = text_from_file("'-deg_wave'"))
l_bkg = Label(lfr_model,text='deg_bkg:', background=bg_frame)
l_bkg.grid(row=4, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_bkg, text = text_from_file("'-deg_bkg'"))

l_rvg = Label(lfr_tpl,text='rvguess:', background=bg_frame)
l_rvg.grid(row=0, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_rvg, text = text_from_file("'-rv_guess'"))
l_overs = Label(lfr_tpl,text='oversampl:', background=bg_frame)
l_overs.grid(row=1, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_overs, text = text_from_file("'-oversampling'"))

l_kaps = Label(lfr_stat,text='kapsig:', background=bg_frame)
l_kaps.grid(row=0, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Help_Box(widget = l_kaps, text = text_from_file("'-kapsig'"))

#Label(fr2,text='stepRV:', background=bg_frame).grid(row=6, column=0, sticky="nw", padx=(xy0,x1), pady=y1)

l_tell = Label(lfr_tell,text='telluric:', background=bg_frame)
l_tell.grid(row=0, column=0, sticky="nw", padx=(xy0,0), pady=y1)
Help_Box(widget = l_tell, text = text_from_file("'-telluric'"))
l_tsig = Label(lfr_tell,text='tsig:', background=bg_frame)
l_tsig.grid(row=1, column=0, sticky="nw", padx=(xy0,0), pady=y1)
Help_Box(widget = l_tsig, text = text_from_file("'-tsig'"))

Label(lfr_tell, text='Optical molecules (TLS, CES, OES, KECK):', background=bg_frame).grid(row=3, column=0, sticky="nw", padx=(xy0,0), pady=y1, columnspan=6)
Label(lfr_tell, text='Near Infrared molecules (CRIRES+):', background=bg_frame).grid(row=5, column=0, sticky="nw", padx=(xy0,0), pady=y1, columnspan=6)

#Label(fr3, text='Plot fitted chunks', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=7, column=0, sticky="nw", padx=(xy0,x1), pady=y1, columnspan=3)

Label(master=win, text='Current command:', background=bg_color).grid(row=2, column = 2, sticky="nw",padx=(0, xy0),pady=10)

###### MAIN ######
FTS_default()

win.mainloop()
