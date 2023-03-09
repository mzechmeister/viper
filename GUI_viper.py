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

###### SETTING ######

font_type = 'Arial'
font_size = 12
bg_frame = '#f1f1f1'    # bg color big frame
bg_color = '#e6e1e1'    # bg color small frames

win_width = 840     # width of GUI window
win_high = 630      # height of GUI window

win = Tk()
win.title('Gui VIPER')
win.geometry("{}x{}".format(win_width, win_high))
win.configure(background=bg_color)


###### POSITIONS ######

x0 = 40     # positions start left
lh = 30     # line height
xs = 100
xss = 80
sp = 210
sps = 170
sl = 7

###### FRAMES #######

# options: solid, sunken, ridge, groove, raised
# groove, ridge need bd=2 or more
Frame(master=win, height=138, width=win_width-40, bg=bg_frame, bd=2, relief='groove').place(x=20, y=20)

Frame(master=win, height=378, width=380, bg=bg_frame, bd=2, relief='groove').place(x=20, y=50+112)
Frame(master=win, height=283, width=414, bg=bg_frame, bd=2, relief='groove').place(x=405, y=50+112)
Frame(master=win, height=115, width=414, bg=bg_frame, bd=2, relief='groove').place(x=405, y=50+112+288)

###### BUTTONS ######

ttk.Style().configure("TButton", padding=2,   background="#fdfdfd", font=(font_type,font_size,'bold'), borderwidth =2)

b_exit = ttk.Button(master=win, text='EXIT', command = bt_exit)
b_exit.place(x=win_width-82,y=win_high-50, width=60)

b_go = ttk.Button(master=win, text='Start', command = bt_start)
b_go.place(x=win_width-162, y=win_high-50, width=60)

b_tpl = Button(master=win, text='Search data file', command = bt_data, background="#fdfdfd")
b_tpl.place(x=win_width-580, y=lh, height=25)

b_tpl = Button(master=win,text='Search tpl file', command = bt_tpl, background="#fdfdfd")
b_tpl.place(x=win_width-168, y=lh, height=25)

###### ENTRIES ######

e_dat = Entry(master=win)
e_dat.insert(0, 'data/TLS/hd189733/*')
e_dat.place(x=x0, y=2*lh, width=350)

e_tpl = Entry(master=win)
e_tpl.insert(0, 'data/TLS/Deconv/HARPS*fits')
e_tpl.place(x=x0+400, y=2*lh, width=350)

e_flag = Entry(master=win)
e_flag.insert(0, '')
e_flag.place(x=x0+xs, y=4*lh, width=340)

e_targ = Entry(master=win)
e_targ.place(x=x0+sp+xss, y=3*lh, width=150)

e_tag = Entry(master=win)
e_tag.place(x=x0+2*sp+xss, y=3*lh, width=240)

e_nset = Entry(master=win)
e_nset.insert(0, ':1')
e_nset.place(x=x0+xss, y=(sl+0)*lh, width=100)

e_oset = Entry(master=win)
e_oset.insert(0, '20')
e_oset.place(x=x0+xss, y=(sl+1)*lh, width=100)

e_ch = Entry(master=win)
e_ch.insert(0, '1')
e_ch.place(x=x0+xss, y=(sl+2)*lh, width=60)

e_deg_norm = Entry(master=win)
e_deg_norm.insert(0, '3')
e_deg_norm.place(x=x0+xss+1*sp, y=(sl+4)*lh, width=50)

e_deg_wave = Entry(master=win)
e_deg_wave.insert(0, '3')
e_deg_wave.place(x=x0+xss+1*sp, y=(sl+5)*lh, width=50)

e_deg_bkg = Entry(master=win)
e_deg_bkg.insert(0, '1')
e_deg_bkg.place(x=x0+xss+1*sp, y=(sl+6)*lh, width=50)

e_tsig = Entry(master=win)
e_tsig.insert(0, '1')
e_tsig.place(x=x0+xss+0*sp, y=(sl+7)*lh, width=50)

e_iphs = Entry(master=win)
e_iphs.insert(0, '50')
e_iphs.place(x=x0+xss+sp, y=(sl+1)*lh, width=50)

e_vg = Entry(master=win)
e_vg.insert(0, '1')
e_vg.place(x=x0+xss+sp, y=(sl+2)*lh, width=50)

e_vcut = Entry(master=win)
e_vcut.insert(0, '100')
e_vcut.place(x=x0+xss+sp, y=(sl+3)*lh, width=50)

e_kapsig = Entry(master=win)
e_kapsig.insert(0, '4.5')
e_kapsig.place(x=x0+xss+0*sp, y=(sl+3)*lh, width=60)

e_overs = Entry(master=win)
e_overs.insert(0, '1')
e_overs.place(x=x0+xss+0*sp, y=(sl+4)*lh, width=60)

e_molec = Entry(master=win)
e_molec.insert(0, 'all')
e_molec.place(x=x0+380, y=(sl+10.5)*lh, width=250)

###### COMBOBOXES ######

combo_inst = ttk.Combobox(master=win, values=['TLS', 'CRIRES','cplCRIRES', 'CES', 'KECK', 'UVES', 'OES'])
combo_inst.set('TLS')
combo_inst.place(x=x0+xs, y=3*lh, width=100)

combo_ip = ttk.Combobox(master=win, values=['g', 'ag', 'agr', 'bg', 'mcg', 'mg', 'sg', 'bnd'])
combo_ip.set('g')
combo_ip.place(x=x0+xss+sp, y=(sl+0)*lh, width=50)

combo_stepRV = ttk.Combobox(master=win, values=['', 'a', 'm'])
combo_stepRV.set('')
combo_stepRV.place(x=x0+xss+0*sp, y=(sl+5)*lh, width=60)

combo_tell = ttk.Combobox(master=win, values=['', 'add', 'add2', 'mask', 'sig'])
combo_tell.set('')
combo_tell.place(x=x0+xss+0*sp, y=(sl+6)*lh, width=60)

###### CHECKBOXES ######

ttk.Style().configure("TCheckbutton", background=bg_frame, bd=0, highlightthickness=0)

cb_flagfile = IntVar()
ttk.Checkbutton(master=win, text=" ", variable=cb_flagfile).place(x=x0, y=4*lh)

cb_wgt = IntVar()
ttk.Checkbutton(master=win, text="     weighted error", variable=cb_wgt).place(x=x0+sp+0*sps, y=(sl+7)*lh)

cb_nocell = IntVar()
ttk.Checkbutton(master=win, text="     no cell", variable=cb_nocell).place(x=x0+sp+0*sps, y=(sl+9)*lh)
cb_createtpl = IntVar()
ttk.Checkbutton(master=win, text="     create tpl", variable=cb_createtpl).place(x=x0, y=(sl+9)*lh)
cb_tellshift = IntVar()
ttk.Checkbutton(master=win, text="     tell shift", variable=cb_tellshift).place(x=x0, y=(sl+10)*lh)

cb_demo = [IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar(), IntVar()]
ttk.Checkbutton(master=win, text="     raw data", variable=cb_demo[0]).place(x=x0+sp+sps, y=(sl+0)*lh)
ttk.Checkbutton(master=win, text="     plot IP", variable=cb_demo[1]).place(x=x0+sp+sps, y=(sl+1)*lh)
ttk.Checkbutton(master=win, text="     stellar tpl", variable=cb_demo[2]).place(x=x0+sp+sps, y=(sl+2)*lh)
ttk.Checkbutton(master=win, text="     forward model", variable=cb_demo[3]).place(x=x0+sp+sps, y=(sl+3)*lh)
ttk.Checkbutton(master=win, text="     fit continuum", variable=cb_demo[4]).place(x=x0+sp+2*sps, y=(sl+0)*lh)
ttk.Checkbutton(master=win, text="     wavelength solution", variable=cb_demo[5]).place(x=x0+sp+2*sps, y=(sl+1)*lh)
ttk.Checkbutton(master=win, text="     fit vguess", variable=cb_demo[6]).place(x=x0+sp+2*sps, y=(sl+2)*lh)
cb_lookguess = IntVar()
ttk.Checkbutton(master=win, text="     lookguess", variable=cb_lookguess).place(x=x0+sp+2*sps, y=(sl+3)*lh)
cb_lookpar = IntVar()
ttk.Checkbutton(master=win, text="     lookpar", variable=cb_lookpar).place(x=x0+sp+2*sps, y=(sl+4)*lh)
cb_lookres = IntVar()
ttk.Checkbutton(master=win, text="     lookres", variable=cb_lookres).place(x=x0+sp+1*sps, y=(sl+4)*lh)
cb_infoprec = IntVar()
ttk.Checkbutton(master=win, text="     infoprec", variable=cb_infoprec).place(x=x0+sp+1*sps, y=(sl+5)*lh)
cb_lookfast = IntVar()
ttk.Checkbutton(master=win, text="     lookfast", variable=cb_lookfast).place(x=x0+sp+1*sps, y=(sl+7)*lh)
cb_lookfast.set(1)
cb_look = IntVar()
ttk.Checkbutton(master=win, text="     look", variable=cb_look).place(x=x0+sp+2*sps, y=(sl+7)*lh)


#cb_atmmask = IntVar()
#ttk.Checkbutton(master=win, text="     atmmask", variable=cb_atmmask).place(x=x0+0*sp, y=(sl+5)*lh)
#cb_atmmod = IntVar()
#ttk.Checkbutton(master=win, text="     atmmod", variable=cb_atmmod).place(x=x0+0*sp, y=(sl+6)*lh)


###### LABELS ######

Label(master=win, text='template file', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0+400,y=lh)
Label(master=win, text='data files', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0, y=lh)
Label(master=win, text='spectrograph:', background=bg_frame).place(x=x0, y=3*lh)
Label(master=win, text='flag file:', background=bg_frame).place(x=x0+20, y=4*lh)
Label(master=win, text='targ:', background=bg_frame).place(x=x0+sp+30, y=3*lh)
Label(master=win, text='tag:', background=bg_frame).place(x=x0+2*sp+30, y=3*lh)

Label(master=win, text='Options data reduction', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0, y=(sl-1.2)*lh)
Label(master=win, text='Options plotting data', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0+sp+sps, y=(sl-1.2)*lh)
Label(master=win, text='Advanced', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0+0*sp, y=(sl+8)*lh)
Label(master=win, text='Molecular specifies (for telluric add)', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0+sp+sps, y=(sl+8.5)*lh)
Label(master=win, text='Optical: H2O O2; NIR: H2O CH4 N2O CO2 CO', background=bg_frame).place(x=x0+sp+sps, y=(sl+9.5)*lh)

Label(master=win, text='nset:', background=bg_frame).place(x=x0, y=(sl+0)*lh)
Label(master=win, text='oset:', background=bg_frame).place(x=x0, y=(sl+1)*lh)
Label(master=win, text='chunks:', background=bg_frame).place(x=x0, y=(sl+2)*lh)

Label(master=win,text='IP:', background=bg_frame).place(x=x0+sp, y=(sl+0)*lh)
Label(master=win,text='iphs:', background=bg_frame).place(x=x0+sp, y=(sl+1)*lh)
Label(master=win,text='rvguess:', background=bg_frame).place(x=x0+sp, y=(sl+2)*lh)
Label(master=win,text='vcut:', background=bg_frame).place(x=x0+1*sp, y=(sl+3)*lh)

Label(master=win,text='deg_norm:', background=bg_frame).place(x=x0+1*sp, y=(sl+4)*lh)
Label(master=win,text='deg_wave:', background=bg_frame).place(x=x0+1*sp, y=(sl+5)*lh)
Label(master=win,text='deg_bkg:', background=bg_frame).place(x=x0+1*sp, y=(sl+6)*lh)

Label(master=win,text='kapsig:', background=bg_frame).place(x=x0+0*sp, y=(sl+3)*lh)
Label(master=win,text='oversampl:', background=bg_frame).place(x=x0+0*sp, y=(sl+4)*lh)
Label(master=win,text='stepRV:', background=bg_frame).place(x=x0+0*sp, y=(sl+5)*lh)

Label(master=win,text='telluric:', background=bg_frame).place(x=x0+0*sp, y=(sl+6)*lh)
Label(master=win,text='tsig:', background=bg_frame).place(x=x0+0*sp, y=(sl+7)*lh)

Label(master=win, text='Plot fitted chunks', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0+sp+sps, y=(sl+6)*lh)

###### MAIN ######

win.mainloop()
