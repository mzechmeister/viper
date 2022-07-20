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

gplot.colors('classic')
gplot2 = Gplot()
g = Gplot()


###### FUNCIONS BUTTONS ######

def bt_start():
    print('---Starting viper.py with choosen parameters---')

    ar_demo = np.zeros(len(cb_demo))
    for d, dem in enumerate(cb_demo):
        ar_demo[d] = dem.get()
    x = np.ones(len(cb_demo))*2    
    demo = (np.poly1d(ar_demo[::-1])(x))[0]

    str_arg = e_dat.get()+"' "+e_tpl.get()+" -inst "+combo_inst.get()+" -dega "+e_dega.get()+" -degb "+e_degb.get()+" -degc "+e_degc.get()+" -nset "+str(e_nset.get())+" -oset "+e_oset.get()+" -chunks "+e_ch.get()+" -demo "+str(int(demo))+" -vg "+e_vg.get() +" -ip "+combo_ip.get() +" -iphs "+e_iphs.get() +" -tsig "+ e_tsig.get()  

#+" -atmmask "+str(cb_atmmask.get()) +" -atmmod "+str(cb_atmmod.get())

    if combo_stepRV.get():
        str_arg += " -stepRV " + str(combo_stepRV.get())
    if combo_tell.get():
        str_arg += " -telluric " + str(combo_tell.get())
    if e_kapsig.get():
        str_arg += " -kapsig "+ str(e_kapsig.get())
    if cb_lookpar.get() == 1:
        str_arg += " -lookpar "
    if cb_lookguess.get() == 1:
        str_arg += " -lookguess "
    if cb_lookres.get() == 1:
        str_arg += " -lookres "
    if cb_infoprec.get() == 1:
        str_arg += " -infoprec "
    if e_overs.get():
        str_arg += " -oversampling "+ str(e_overs.get())
    if e_targ.get():
        str_arg += " -targ " + str(e_targ.get())
    if e_tag.get():
        str_arg += " -tag " + str(e_tag.get())

    os.system("python3 viper.py '" + str_arg)

   # doesn't work as subprocess
   # os.system("python3 viper.py 'data/CRIRES/Ross619/210*' data/CRIRES/model/Ross619_19.fits -inst CRIRES -dega 2 -degb 2 -nset :2 -oset 10")
  #  subprocess.Popen(["./viper.py","'data/CRIRES/Ross619/210*'", "data/CRIRES/model/Ross619_19.fits -inst CRIRES -dega 2 -degb 2 -nset :2 -oset 10 -test 1"])
  
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
bg_frame = '#f1f1f1'	# bg color big frame
bg_color = '#e6e1e1'	# bg color small frames

win_width = 840		# width of GUI window
win_high = 410		# height of GUI window

win = Tk()
win.title('Gui VIPER')
win.geometry("{}x{}".format(win_width, win_high))
win.configure(background=bg_color)


###### POSITIONS ######

x0 = 40     # positions start left
lh = 30     # line height
xs = 100
xss = 70
sp = 220
sps = 170
sl = 6

###### FRAMES #######

# options: solid, sunken, ridge, groove, raised
# groove, ridge need bd=2 or more
Frame(master=win,height=108,width=win_width-40,bg=bg_frame,bd=2, relief='groove').place(x=20,y=20)

Frame(master=win,height=258,width=380,bg=bg_frame,bd=2, relief='groove').place(x=20,y=20+112)
Frame(master=win,height=223,width=414,bg=bg_frame,bd=2, relief='groove').place(x=405,y=20+112)


###### BUTTONS ######

ttk.Style().configure("TButton", padding=2,   background="#fdfdfd",font=(font_type,font_size,'bold'),borderwidth =2)

b_exit = ttk.Button(master=win,text='EXIT', command = bt_exit)
b_exit.place(x=win_width-82,y=win_high-50, width=60)

b_go = ttk.Button(master=win,text='Start', command = bt_start)
b_go.place(x=win_width-162,y=win_high-50, width=60)

b_tpl = Button(master=win,text='Search data file', command = bt_data, background="#fdfdfd")
b_tpl.place(x=win_width-580,y=lh, height=25)

b_tpl = Button(master=win,text='Search tpl file', command = bt_tpl, background="#fdfdfd")
b_tpl.place(x=win_width-168,y=lh, height=25)

###### ENTRIES ######

e_dat = Entry(master=win)
e_dat.insert(0,'data/TLS/hd189733/*')
e_dat.place(x=x0,y=2*lh,width=350)

e_tpl = Entry(master=win)
e_tpl.insert(0,'data/TLS/Deconv/HARPS*fits')
e_tpl.place(x=x0+400,y=2*lh,width=350)

e_targ = Entry(master=win)
e_targ.place(x=x0+sp+xss,y=3*lh,width=150)

e_tag = Entry(master=win)
e_tag.place(x=x0+2*sp+xss,y=3*lh,width=240)

e_nset = Entry(master=win)
e_nset.insert(0,':1')
e_nset.place(x=x0+xss,y=(sl+0)*lh,width=100)

e_oset = Entry(master=win)
e_oset.insert(0,'20')
e_oset.place(x=x0+xss,y=(sl+1)*lh,width=100)

e_ch = Entry(master=win)
e_ch.insert(0,'1')
e_ch.place(x=x0+xss,y=(sl+2)*lh,width=60)

e_dega = Entry(master=win)
e_dega.insert(0,'3')
e_dega.place(x=x0+xss+1*sp,y=(sl+3)*lh,width=50)

e_degb = Entry(master=win)
e_degb.insert(0,'3')
e_degb.place(x=x0+xss+1*sp,y=(sl+4)*lh,width=50)

e_degc = Entry(master=win)
e_degc.insert(0,'1')
e_degc.place(x=x0+xss+1*sp,y=(sl+5)*lh,width=50)

e_tsig = Entry(master=win)
e_tsig.insert(0,'1')
e_tsig.place(x=x0+xss+1*sp,y=(sl+6)*lh,width=50)

e_iphs = Entry(master=win)
e_iphs.insert(0,'50')
e_iphs.place(x=x0+xss+sp,y=(sl+1)*lh,width=50)

e_vg = Entry(master=win)
e_vg.insert(0,'1')
e_vg.place(x=x0+xss+sp,y=(sl+2)*lh,width=50)

e_kapsig = Entry(master=win)
e_kapsig.insert(0,'4.5')
e_kapsig.place(x=x0+xss+0*sp,y=(sl+3)*lh,width=60)

e_overs = Entry(master=win)
e_overs.insert(0,'1')
e_overs.place(x=x0+xss+0*sp,y=(sl+4)*lh,width=60)

###### COMBOBOXES ######

combo_inst = ttk.Combobox(master=win,values=['TLS','CRIRES','CES', 'KECK', 'UVES','OES'])
combo_inst.set('TLS')
combo_inst.place(x=x0+xs,y=3*lh,width=100)

combo_ip = ttk.Combobox(master=win,values=['g', 'ag', 'sg', 'mg', 'bnd'])
combo_ip.set('g')
combo_ip.place(x=x0+xss+sp,y=(sl+0)*lh,width=50)

combo_stepRV = ttk.Combobox(master=win,values=['', 'a', 'm'])
combo_stepRV.set('')
combo_stepRV.place(x=x0+xss+0*sp,y=(sl+5)*lh,width=60)

combo_tell = ttk.Combobox(master=win,values=['','add', 'mask', 'sig'])
combo_tell.set('')
combo_tell.place(x=x0+xss+0*sp,y=(sl+6)*lh,width=60)

###### CHECKBOXES ######

ttk.Style().configure("TCheckbutton",background =bg_frame,bd=0, highlightthickness=0)

cb_demo = [IntVar(),IntVar(),IntVar(),IntVar(),IntVar(),IntVar(),IntVar()]

ttk.Checkbutton(master=win, text="     raw data", variable=cb_demo[0]).place(x=x0+sp+sps,y=(sl+0)*lh)
ttk.Checkbutton(master=win, text="     plot IP", variable=cb_demo[1]).place(x=x0+sp+sps,y=(sl+1)*lh)
ttk.Checkbutton(master=win, text="     stellar tpl", variable=cb_demo[2]).place(x=x0+sp+sps,y=(sl+2)*lh)
ttk.Checkbutton(master=win, text="     forward model", variable=cb_demo[3]).place(x=x0+sp+sps,y=(sl+3)*lh)
ttk.Checkbutton(master=win, text="     fit continuum", variable=cb_demo[4]).place(x=x0+sp+2*sps,y=(sl+0)*lh)
ttk.Checkbutton(master=win, text="     wavelength solution", variable=cb_demo[5]).place(x=x0+sp+2*sps,y=(sl+1)*lh)
ttk.Checkbutton(master=win, text="     fit vguess", variable=cb_demo[6]).place(x=x0+sp+2*sps,y=(sl+2)*lh)
cb_lookguess = IntVar()
ttk.Checkbutton(master=win, text="     lookguess", variable=cb_lookguess).place(x=x0+sp+2*sps,y=(sl+3)*lh)
cb_lookpar = IntVar()
ttk.Checkbutton(master=win, text="     lookpar", variable=cb_lookpar).place(x=x0+sp+2*sps,y=(sl+4)*lh)
cb_lookres = IntVar()
ttk.Checkbutton(master=win, text="     lookres", variable=cb_lookres).place(x=x0+sp+1*sps,y=(sl+4)*lh)
cb_infoprec = IntVar()
ttk.Checkbutton(master=win, text="     infoprec", variable=cb_infoprec).place(x=x0+sp+1*sps,y=(sl+5)*lh)

#cb_atmmask = IntVar()
#ttk.Checkbutton(master=win, text="     atmmask", variable=cb_atmmask).place(x=x0+0*sp,y=(sl+5)*lh)
#cb_atmmod = IntVar()
#ttk.Checkbutton(master=win, text="     atmmod", variable=cb_atmmod).place(x=x0+0*sp,y=(sl+6)*lh)


###### LABELS ######

(Label(master=win,text='template file', font=(font_type,font_size,'bold'),background=bg_frame)).place(x=x0+400,y=lh)
(Label(master=win,text='data files', font=(font_type,font_size,'bold'),background=bg_frame)).place(x=x0,y=lh)
(Label(master=win,text='telescope:', background=bg_frame)).place(x=x0,y=3*lh)
(Label(master=win,text='targ:', background=bg_frame)).place(x=x0+sp+30,y=3*lh)
(Label(master=win,text='tag:', background=bg_frame)).place(x=x0+2*sp+30,y=3*lh)

(Label(master=win,text='Options data reduction', font=(font_type,font_size,'bold'),background=bg_frame)).place(x=x0,y=5*lh-5)
(Label(master=win,text='Options plotting data', font=(font_type,font_size,'bold'), background=bg_frame)).place(x=x0+sp+sps,y=5*lh-5)

(Label(master=win,text='nset:', background=bg_frame)).place(x=x0,y=(sl+0)*lh)
(Label(master=win,text='oset:', background=bg_frame)).place(x=x0,y=(sl+1)*lh)
(Label(master=win,text='chunks:', background=bg_frame)).place(x=x0,y=(sl+2)*lh)

(Label(master=win,text='dega:', background=bg_frame)).place(x=x0+1*sp,y=(sl+3)*lh)
(Label(master=win,text='degb:', background=bg_frame)).place(x=x0+1*sp,y=(sl+4)*lh)
(Label(master=win,text='degc:', background=bg_frame)).place(x=x0+1*sp,y=(sl+5)*lh)

(Label(master=win,text='IP:', background=bg_frame)).place(x=x0+sp,y=(sl+0)*lh)
(Label(master=win,text='iphs:', background=bg_frame)).place(x=x0+sp,y=(sl+1)*lh)
(Label(master=win,text='vguess:', background=bg_frame)).place(x=x0+sp,y=(sl+2)*lh)

(Label(master=win,text='kapsig:', background=bg_frame)).place(x=x0+0*sp,y=(sl+3)*lh)
(Label(master=win,text='oversampl:', background=bg_frame)).place(x=x0+0*sp,y=(sl+4)*lh)
(Label(master=win,text='stepRV:', background=bg_frame)).place(x=x0+0*sp,y=(sl+5)*lh)

(Label(master=win,text='telluric:', background=bg_frame)).place(x=x0+0*sp,y=(sl+6)*lh)
(Label(master=win,text='tsig:', background=bg_frame)).place(x=x0+1*sp,y=(sl+6)*lh)


###### MAIN ######

win.mainloop()
