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

    str_arg = e_dat.get()+"' "+e_tpl.get()+" -inst "+combo_inst.get()+" -deg_norm "+e_deg_norm.get()+" -deg_wave "+e_deg_wave.get()+" -deg_bkg "+e_deg_bkg.get()+" -nset "+str(e_nset.get())+" -oset "+e_oset.get()+" -chunks "+e_ch.get()+" -demo "+str(int(demo))+" -rv_guess "+e_vg.get() +" -ip "+combo_ip.get() +" -iphs "+e_iphs.get() +" -tsig "+ e_tsig.get() +" -vcut "+e_vcut.get() +" -molec "+e_molec.get() 

#+" -atmmask "+str(cb_atmmask.get()) +" -atmmod "+str(cb_atmmod.get())

    if combo_stepRV.get():
        str_arg += " -stepRV " + str(combo_stepRV.get())
    if combo_tell.get():
        str_arg += " -telluric " + str(combo_tell.get())
    if e_kapsig.get():
        str_arg += " -kapsig "+ str(e_kapsig.get())
    if cb_lookpar.get():
        str_arg += " -lookpar "
    if cb_look.get():
        str_arg += " -look "
    if cb_lookfast.get():
        str_arg += " -lookfast "
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
    if cb_nocell.get():
        str_arg += " -nocell "
    if cb_createtpl.get():
        str_arg += " -createtpl "
    if cb_tellshift.get():
        str_arg += " -tellshift "
    if cb_wgt.get():
        str_arg += " -wgt "

    print(str_arg)

    os.system("python3 "+viperdir+"viper.py '" + str_arg)

   # doesn't work as subprocess
   # os.system("python3 viper.py 'data/CRIRES/Ross619/210*' data/CRIRES/model/Ross619_19.fits -inst CRIRES -deg_norm 2 -deg_wave 2 -nset :2 -oset 10")
  #  subprocess.Popen(["./viper.py","'data/CRIRES/Ross619/210*'", "data/CRIRES/model/Ross619_19.fits -inst CRIRES -deg_norm 2 -deg_wave 2 -nset :2 -oset 10 -test 1"])

    print('---Finished viper.py---')

def bt_data():
    data_file = askopenfilename(master=win)
    if data_file:
        e_dat.delete(0, END)
        e_dat.insert(0, data_file)

def bt_tpl():
    tpl_file = askopenfilename(master=win)
    if tpl_file:
        e_tpl.delete(0, END)
        e_tpl.insert(0, tpl_file)

def bt_exit():
    exit()

def new(event):
    win_width = win.winfo_width()
    win_height = event.height
    fr1.config(width=win_width-40)
    fr2.config(width=(win_width-46)/2)
    fr3.config(width=(win_width-46)/2)
    fr4.config(width=(win_width-46)/2)

###### SETTING ######

font_type = 'Arial'
font_size = 12
bg_frame = '#f1f1f1'    # bg color big frame
bg_color = '#e6e1e1'    # bg color small frames

win_width = 940     # width of GUI window	# 840
win_high = 690      # height of GUI window	# 630
xy0 = 20
y1 = 5
x1 = 10

win = Tk()
win.title('Gui VIPER')
win.geometry("{}x{}".format(win_width, win_high))
win.configure(background=bg_color)


###### FRAMES #######
win.bind("<Configure>", new)
win.columnconfigure(0, weight=1)
win.columnconfigure(3, weight=1)
#win.rowconfigure(0, weight=1)

# options: solid, sunken, ridge, groove, raised
# groove, ridge need bd=2 or more
fr1 = Frame(win, height=148, width=win_width-40, bg=bg_frame, bd=2, relief='groove')
fr1.grid(row=0, column = 1, sticky="nw",padx=20,pady=(xy0,6), columnspan=2)
fr1.grid_propagate(False)
fr1.grid_columnconfigure(2, weight=1)
fr1.grid_columnconfigure(4, weight=1)
fr1.grid_rowconfigure(0, weight=1)
fr1.grid_rowconfigure(5, weight=1)

fr2 = Frame(master=win, height=428, width=(win_width-46)/2, bg=bg_frame, bd=2, relief='groove')
fr2.grid(row=1, column = 1, sticky="nw",padx=(xy0,6),pady=0, rowspan=2)
fr2.grid_propagate(False)
fr2.grid_columnconfigure(1, weight=1)
fr2.grid_columnconfigure(3, weight=1)
#fr2.grid_rowconfigure(0, weight=1)

fr3 = Frame(master=win, height=303, width=(win_width-46)/2, bg=bg_frame, bd=2, relief='groove')
fr3.grid(row=1, column = 2, sticky="nw",padx=(0, xy0),pady=0)
fr3.grid_propagate(False)
fr3.grid_columnconfigure(0, weight=1)
fr3.grid_columnconfigure(1, weight=1)

fr4 = Frame(master=win, height=135, width=(win_width-46)/2, bg=bg_frame, bd=2, relief='groove')
fr4.grid(row=2, column = 2, sticky="nw",padx=(0, xy0),pady=(6,0))
fr4.grid_propagate(False)

###### BUTTONS ######

ttk.Style().configure("TButton", padding=2,   background="#fdfdfd", font=(font_type,font_size,'bold'), borderwidth =2)

b_exit = ttk.Button(master=win, text='EXIT', command = bt_exit)
b_exit.grid(row=3, column=2, sticky="ne", padx=(0,xy0), pady = xy0)

b_go = ttk.Button(master=win, text='Start', command = bt_start)
b_go.grid(row=3, column=2, sticky="ne", padx=(0,130), pady = xy0)

b_dat = Button(fr1, text='Search data file', command = bt_data, background="#fdfdfd", width=15)
b_dat.grid(row=1, column=7, sticky="nw", padx=xy0)

b_tpl = Button(fr1,text='Search tpl file', command = bt_tpl, background="#fdfdfd", width=15)
b_tpl.grid(row=2, column=7, sticky="nw", padx=xy0)

###### ENTRIES ######

e_dat = Entry(fr1, width = 200)
e_dat.insert(0, 'data/TLS/hd189733/*')
e_dat.grid(row=1, column=1, sticky="nw", padx=x1, pady=y1, columnspan=6)

e_tpl = Entry(fr1, width = 200)
e_tpl.insert(0, 'data/TLS/Deconv/HARPS*fits')
e_tpl.grid(row=2, column=1, sticky="nw", padx=x1, pady=y1, columnspan=6)

e_flag = Entry(fr1, width = 200)
e_flag.insert(0, '')
e_flag.grid(row=3, column=1, sticky="nw", padx=x1, pady=y1, columnspan=6)

e_targ = Entry(fr1, width=70)
e_targ.grid(row=4, column=4, sticky="nw", padx=x1, pady=y1)

e_tag = Entry(fr1, width=35)
e_tag.grid(row=4, column=6, sticky="nw", padx=xy0, pady=y1, columnspan=2)

e_nset = Entry(fr2)#, width=15)
e_nset.insert(0, ':1')
e_nset.grid(row=1, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

e_oset = Entry(fr2)
e_oset.insert(0, '20')
e_oset.grid(row=2, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

e_ch = Entry(fr2)
e_ch.insert(0, '1')
e_ch.grid(row=3, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

e_deg_norm = Entry(fr2)
e_deg_norm.insert(0, '3')
e_deg_norm.grid(row=5, column=3, sticky="nw", padx=(x1,xy0), pady=(y1,0))

e_deg_wave = Entry(fr2)
e_deg_wave.insert(0, '3')
e_deg_wave.grid(row=6, column=3, sticky="nw", padx=(x1,xy0), pady=(y1,0))

e_deg_bkg = Entry(fr2)
e_deg_bkg.insert(0, '1')
e_deg_bkg.grid(row=7, column=3, sticky="nw", padx=(x1,xy0), pady=(y1,0))

e_tsig = Entry(fr2)
e_tsig.insert(0, '1')
e_tsig.grid(row=8, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

e_iphs = Entry(fr2)
e_iphs.insert(0, '50')
e_iphs.grid(row=2, column=3, sticky="nw", padx=(x1,xy0), pady=(y1,0))

e_vg = Entry(fr2)
e_vg.insert(0, '1')
e_vg.grid(row=3, column=3, sticky="nw", padx=(x1,xy0), pady=(y1,0))

e_vcut = Entry(fr2)
e_vcut.insert(0, '100')
e_vcut.grid(row=4, column=3, sticky="nw", padx=(x1,xy0), pady=(y1,0))

e_kapsig = Entry(fr2)
e_kapsig.insert(0, '4.5')
e_kapsig.grid(row=4, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

e_overs = Entry(fr2)
e_overs.insert(0, '1')
e_overs.grid(row=5, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

e_molec = Entry(fr4)
e_molec.insert(0, 'all')
e_molec.grid(row=2, column=0, sticky="nw", padx=(xy0,0), pady=y1, columnspan=3)

###### COMBOBOXES ######

combo_inst = ttk.Combobox(fr1, values=['TLS', 'CRIRES','cplCRIRES', 'CES', 'KECK', 'UVES', 'OES'], width=10)
combo_inst.set('TLS')
combo_inst.grid(row=4, column=1, sticky="nw", padx=x1, pady=y1)

combo_ip = ttk.Combobox(fr2, values=['g', 'ag', 'agr', 'bg', 'mcg', 'mg', 'sg', 'bnd'])
combo_ip.set('g')
combo_ip.grid(row=1, column=3, sticky="nw", padx=(x1,xy0), pady=y1)

combo_stepRV = ttk.Combobox(fr2, values=['', 'a', 'm'])
combo_stepRV.set('')
combo_stepRV.grid(row=6, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

combo_tell = ttk.Combobox(fr2, values=['', 'add', 'add2', 'mask', 'sig'])
combo_tell.set('')
combo_tell.grid(row=7, column=1, sticky="nw", padx=(x1,xy0), pady=y1)

###### CHECKBOXES ######

ttk.Style().configure("TCheckbutton", background=bg_frame, bd=0, highlightthickness=0)

cb_flagfile = IntVar()
ttk.Checkbutton(fr1, text="   flag file:", variable=cb_flagfile).grid(row=3, column=0, sticky="nw", padx=x1, pady=y1)

cb_wgt = IntVar()
ttk.Checkbutton(fr2, text="     weighted error", variable=cb_wgt).grid(row=8, column=2, sticky="nw", padx=(xy0,x1), pady=y1, columnspan=2)

cb_nocell = IntVar()
ttk.Checkbutton(fr2, text="     no cell", variable=cb_nocell).grid(row=x1, column=2, sticky="nw", padx=(xy0,x1), pady=y1)
cb_createtpl = IntVar()
ttk.Checkbutton(fr2, text="     create tpl", variable=cb_createtpl).grid(row=x1, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
cb_tellshift = IntVar()
ttk.Checkbutton(fr2, text="     tell shift", variable=cb_tellshift).grid(row=11, column=0, sticky="nw", padx=(xy0,x1), pady=y1)

cb_demo = [IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar()]
ttk.Checkbutton(fr3, text="     raw data", variable=cb_demo[0]).grid(row=1, column=0, sticky="nw", padx=(xy0,x1), pady=y1)

ttk.Checkbutton(fr3, text="     plot IP", variable=cb_demo[1]).grid(row=2, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(fr3, text="     stellar tpl", variable=cb_demo[2]).grid(row=3, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(fr3, text="     forward model", variable=cb_demo[3]).grid(row=4, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(fr3, text="     fit continuum", variable=cb_demo[4]).grid(row=1, column=1, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(fr3, text="     wavelength solution", variable=cb_demo[5]).grid(row=2, column=1, sticky="nw", padx=(xy0,x1), pady=y1)
ttk.Checkbutton(fr3, text="     fit vguess", variable=cb_demo[6]).grid(row=3, column=1, sticky="nw", padx=(xy0,x1), pady=y1)
cb_lookguess = IntVar()
ttk.Checkbutton(fr3, text="     lookguess", variable=cb_lookguess).grid(row=4, column=1, sticky="nw", padx=(xy0,x1), pady=y1)
cb_lookpar = IntVar()
ttk.Checkbutton(fr3, text="     lookpar", variable=cb_lookpar).grid(row=5, column=1, sticky="nw", padx=(xy0,x1), pady=y1)
cb_lookres = IntVar()
ttk.Checkbutton(fr3, text="     lookres", variable=cb_lookres).grid(row=5, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
cb_infoprec = IntVar()
ttk.Checkbutton(fr3, text="     infoprec", variable=cb_infoprec).grid(row=6, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
cb_lookfast = IntVar()
ttk.Checkbutton(fr3, text="     lookfast", variable=cb_lookfast).grid(row=8, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
cb_lookfast.set(1)
cb_look = IntVar()
ttk.Checkbutton(fr3, text="     look", variable=cb_look).grid(row=8, column=1, sticky="nw", padx=(xy0,x1), pady=y1)


###### LABELS ######


Label(fr1, text='data files', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=1, column=0, sticky="nw", padx=x1, pady=y1)
Label(fr1, text='template file', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=2, column=0, sticky="nw", padx=x1, pady=y1)

Label(fr1, text='spectrograph:', background=bg_frame).grid(row=4, column=0, sticky="nw", padx=x1, pady=y1)

Label(fr1, text='targ:', background=bg_frame).grid(row=4, column=3, sticky="nw", padx=xy0, pady=y1)
Label(fr1, text='tag:', background=bg_frame).grid(row=4, column=5, sticky="nw", padx=xy0, pady=y1)

Label(fr2, text='Options data reduction', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=0, column=0, sticky="nw", padx=(xy0,0), pady=(xy0,y1), columnspan=3)
Label(fr3, text='Options plotting data', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=0, column=0, sticky="nw", padx=(xy0,0), pady=(xy0,y1), columnspan=3)
Label(fr2, text='Advanced', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=9, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Label(fr4, text='Molecular specifies (for telluric add)', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=0, column=0, sticky="nw", padx=(xy0,0), pady=(xy0,y1), columnspan=3)
Label(fr4, text='Optical: H2O O2; NIR: H2O CH4 N2O CO2 CO', background=bg_frame).grid(row=1, column=0, sticky="nw", padx=(xy0,0), pady=y1, columnspan=3)

Label(fr2, text='nset:', background=bg_frame).grid(row=1, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Label(fr2, text='oset:', background=bg_frame).grid(row=2, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Label(fr2, text='chunks:', background=bg_frame).grid(row=3, column=0, sticky="nw", padx=(xy0,x1), pady=y1)

Label(fr2,text='IP:', background=bg_frame).grid(row=1, column=2, sticky="nw", padx=x1, pady=y1)
Label(fr2,text='iphs:', background=bg_frame).grid(row=2, column=2, sticky="nw", padx=x1, pady=y1)
Label(fr2,text='rvguess:', background=bg_frame).grid(row=3, column=2, sticky="nw", padx=x1, pady=y1)
Label(fr2,text='vcut:', background=bg_frame).grid(row=4, column=2, sticky="nw", padx=x1, pady=y1)

Label(fr2,text='deg_norm:', background=bg_frame).grid(row=5, column=2, sticky="nw", padx=x1, pady=y1)
Label(fr2,text='deg_wave:', background=bg_frame).grid(row=6, column=2, sticky="nw", padx=x1, pady=y1)
Label(fr2,text='deg_bkg:', background=bg_frame).grid(row=7, column=2, sticky="nw", padx=x1, pady=y1)

Label(fr2,text='kapsig:', background=bg_frame).grid(row=4, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Label(fr2,text='oversampl:', background=bg_frame).grid(row=5, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Label(fr2,text='stepRV:', background=bg_frame).grid(row=6, column=0, sticky="nw", padx=(xy0,x1), pady=y1)

Label(fr2,text='telluric:', background=bg_frame).grid(row=7, column=0, sticky="nw", padx=(xy0,x1), pady=y1)
Label(fr2,text='tsig:', background=bg_frame).grid(row=8, column=0, sticky="nw", padx=(xy0,x1), pady=y1)

Label(fr3, text='Plot fitted chunks', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=7, column=0, sticky="nw", padx=(xy0,x1), pady=y1, columnspan=3)

###### MAIN ######

win.mainloop()
