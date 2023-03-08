#! /usr/bin/env python3
# -*- coding: iso-8859-1 -*-

# GUI to run vpr.py, showing most important options

import os
from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText

import numpy as np
import sys
import vpr
vpr.pause = print   # mainloop of the gui will pause

###### FUNCIONS BUTTONS ######

options = ['-plot rv', '-plot rvo', '-plot nrvo']
opt = 0

def update_changes(*args, refresh = 0):
    # Update when changing checkboxes
    if refresh:
        refresh_oset(str(refresh))
    if opt < 3:
        call_vpr(options[opt])
    elif opt == 3:
        call_vpr(cmp=True)

def call_vpr(args='-plot rv', cmp=False, res=False):
    print('---Starting vpr.py with choosen parameters---')

    global opt
    if '-save' in str(args):    
        args = options[opt] + ' -save'
    else:
        opt = options.index(str(args))

    if not res:
        str_arg = e_rvo1.get() + ' -avg ' + combo_avg.get()

        if cb_cen.get():
            str_arg += " -cen "
        if cb_ocen.get():
            str_arg += " -ocen "
        if cb_orders1 != []:
            str_arg += " -oset "
            for i,o in enumerate (cb_orders1):
                if o.get():
                    str_arg += str(o_rvo1[i])+","
        # if e_nset.get():
        #     str_arg += " -nset "+str(e_nset.get())
        if e_sort.get():
            str_arg += " -sort "+str(e_sort.get())
        if e_offset.get():
            str_arg += " -offset "+str(e_offset.get())

        if cmp:
            opt = 3
            if e_rvo2.get():
               str_arg += " -cmp "+e_rvo2.get()
            if cb_orders2 != []:
               str_arg += " -cmposet "
               for i,o in enumerate (cb_orders2):
                    if o.get():
                        str_arg += str(o_rvo2[i])+","
            if cb_cmpocen.get():
                str_arg += " -cmpocen"

    else:
        str_arg = " -res "+str(e_dir.get())

        if e_oset_r.get():
            str_arg += " -oset "+e_oset_r.get()
        if e_nset_r.get():
            str_arg += " -nset "+e_nset_r.get()

    if args:
        str_arg += " " + args

    e_run.delete('0.0', END)
    e_run.insert(INSERT,"python3 vpr.py "+str_arg)

    print("option string:", str_arg)
    vpr.run(str_arg.split())

    print('---Finished vpr.py---')

def bt_rvo1():
    file1 = askopenfilename()
    if file1:
        filename1.set(file1)
        refresh_oset('1')
        update_changes()

def bt_rvo2():
    file2 = askopenfilename()
    if file2:
        filename2.set(file2)
        refresh_oset('2')
        update_changes()

def bt_swap():
    f1 = filename1.get()
    f2 = filename2.get()
    filename1.set(f2)
    filename2.set(f1)
    refresh_oset('12')
    update_changes()


cbo1, cbo2 = [], []			# oset checkboxes
cb_orders1, cb_orders2 = [], []		# Var() of  oset checkboxes 
o_rvo1, o_rvo2 = [], []			# orders in RV files

def refresh_oset(num):
    # refresh oset checkboxes (when using new RV file)
    if '1' in num:
        if e_rvo1.get():
            global cb_orders1, cbo1, o_rvo1 
            o_rvo1 = get_orders(e_rvo1.get())
            cb_orders1, cbo1 = create_cb(o_rvo1, cbo1, cb_orders1, 0)
    if '2' in num:
        if e_rvo2.get():
            global cb_orders2, cbo2, o_rvo2 
            o_rvo2 = get_orders(e_rvo2.get())
            cb_orders2, cbo2 = create_cb(o_rvo2, cbo2, cb_orders2, 1)

def set_oset(cb_orders, value):
    # set all oset variables to one value
    if value == 0:
        global opt
        opt = 10
    for o in cb_orders:
        o.set(value)

def get_orders(file):
    # get all available orders from RV file
    Afull = np.genfromtxt(file, dtype=None, names=True,
            deletechars='',   # to keep the dash for chunks
            encoding=None).view(np.recarray)
    colnames = Afull.dtype.names

    onames_all = np.array([col for col in colnames if col.startswith('rv')], dtype='O')
    orders_all, chunks_all = [*np.array([[*map(int, oname[2:].split('-'))] for oname in onames_all]).T, None][0:2]

    return orders_all

def create_cb(orders_all, cbo, cb_orders, pos):
    # creates checkboxes for oset
    if cbo != []:
        # delete old checkboxes
        for i in cbo:
            i.destroy()
        cb_orders = []
        cbo = []    

    for i, o in enumerate(orders_all):
        cb_orders.append(IntVar())
        yi, xi =  divmod(i, 10)
        c = ttk.Checkbutton(tab_rv, text="  "+str(o),  variable=cb_orders[i])
        c.place(x=x0+75*xi, y=(6+yi+(pos*4))*lh*0.8)
        cbo.append(c)
        cb_orders[i].set(1)
        cb_orders[i].trace("w", update_changes)

    return cb_orders, cbo

def bt_exit():
    exit()


###### SETTING ######

font_type = 'Arial'
font_size = 12
bg_frame = '#f1f1f1'    # bg color big frame
bg_color = '#e6e1e1'    # bg color small frames

win_width = 840     # width of GUI window
win_high = 655      # height of GUI window

win = Tk()
win.title('Gui VPR')
win.geometry("{}x{}".format(win_width, win_high))
win.configure(background=ttk.Style().lookup('TFrame', 'background'))   # same bg-color for win and Frame

style = ttk.Style(win)
style.configure('lefttab.TNotebook', tabposition='nw')
style.configure('TNotebook.Tab', font=(font_type, font_size+2, 'bold'))

tabControl = ttk.Notebook(win, style='lefttab.TNotebook')

tab_rv = ttk.Frame(tabControl)
tab_res = ttk.Frame(tabControl)

tabControl.add(tab_rv, text ='RV plots')
tabControl.add(tab_res, text ='Residuals')
tabControl.pack(expand=1, fill="both", side='top')


###### POSITIONS ######

x0 = 40     # positions start left
lh = 30     # line height
xs = 100
xss = 80
sp = 270
sps = 170
sl = 7
fr_high = win_high-170
width_bt = 110

###### FRAMES #######

# options: solid, sunken, ridge, groove, raised
# groove, ridge need bd=2 or more
frm_rv = Frame(tab_rv, bg=bg_frame, bd=2, relief='groove')
frm_rv.pack(padx=5, pady=2, expand=True, fill='both', side='left')
Frame(tab_res, height=fr_high, width=win_width-40, bg=bg_frame, bd=2, relief='groove').pack(padx=5, pady=2, expand=True, fill='both')

#Frame(tab_rv, height=185, width=390, bg=bg_frame, bd=2, relief='groove').place(x=20, y=50+112)
#Frame(tab_rv, height=185, width=404, bg=bg_frame, bd=2, relief='groove').place(x=415, y=50+112)

file_bar = Frame(frm_rv, bg=bg_frame)
file_bar.pack(padx=20, pady=(lh*2, 0), fill='x', side='top')
plot_bar = Frame(frm_rv, bg=bg_frame)
plot_bar.pack(padx=20, pady=(lh*12, 0), fill='x', side='top')

###### BUTTONS ######

ttk.Style().configure("TButton", padding=2,   background="#fdfdfd", font=(font_type,font_size,'bold'), borderwidth =2)
ttk.Style().configure("B2.TButton", padding=1,   background="#fdfdfd", font=(font_type,font_size-1,''), borderwidth =1)

b_exit = ttk.Button(master=win, text='EXIT', command = bt_exit)
b_exit.pack(side='right', padx=30)

b_rvo1 = Button(tab_rv, text='Search data file', command = bt_rvo1, background="#fdfdfd")
b_rvo1.place(x=win_width-590, y=lh, height=25)

b_swap = Button(file_bar, text='swap', command = bt_swap, background="#fdfdfd")

b_rvo2 = Button(tab_rv,text='Search data file', command = bt_rvo2, background="#fdfdfd")
b_rvo2.place(x=win_width-169, y=lh, height=25)

conf = {'side': 'left', 'padx': 10}

b_rvbjd = ttk.Button(plot_bar, text='Plot BJD-RV', command = call_vpr)
b_rvbjd.pack(**conf)

b_rvo = ttk.Button(plot_bar, text='Plot o-rv', command = lambda: call_vpr('-plot rvo'))
b_rvo.pack(**conf)

b_nrvo = ttk.Button(plot_bar, text='Plot no-rv', command = lambda: call_vpr('-plot nrvo'))
b_nrvo.pack(**conf)

b_cmp = ttk.Button(plot_bar, text='Compare', command = lambda: call_vpr(cmp=True))
b_cmp.pack(**{**conf, 'side': 'right'})

b_save = ttk.Button(tab_rv, text='Save', command = lambda: call_vpr('-save'))
b_save.place(x=x0+300, y=(13.5)*lh-5, width=width_bt/2)

b_res = ttk.Button(tab_res, text='Plot res', command = lambda: call_vpr(res=True))
b_res.place(x=win_width-(width_bt+20)*1-20, y=fr_high-30, width=width_bt)


b_oset1_a = ttk.Button(tab_rv, text='select all', style = 'B2.TButton', command = lambda: set_oset(cb_orders1, 1))
b_oset1_a.place(x=x0+100, y=4*lh, width=width_bt*0.8, height=22)
b_oset1_n = ttk.Button(tab_rv, text='select none', style = 'B2.TButton', command = lambda: set_oset(cb_orders1, 0))
b_oset1_n.place(x=x0+200, y=4*lh, width=width_bt*0.8, height=22)

b_oset2_a = ttk.Button(tab_rv, text='select all', style = 'B2.TButton', command = lambda: set_oset(cb_orders2, 1))
b_oset2_a.place(x=x0+100, y=7*lh, width=width_bt*0.8, height=22)
b_oset2_n = ttk.Button(tab_rv, text='select none', style = 'B2.TButton', command = lambda: set_oset(cb_orders2, 0))
b_oset2_n.place(x=x0+200, y=7*lh, width=width_bt*0.8, height=22)


###### ENTRIES ######

filename1 = StringVar()
filename2 = StringVar()

e_run = ScrolledText(win, background='#f0f0f0')

e_rvo1 = Entry(file_bar, textvariable=filename1, width=40)
e_rvo1.insert(0, 'tmp.rvo.dat')
e_rvo1.bind("<Return>", (lambda event: update_changes(refresh=1)))
e_rvo1.pack(side='left', padx=20, pady=3, expand=1, fill='x')

b_swap.pack(side='left')

e_rvo2 = Entry(file_bar, textvariable=filename2)
e_rvo2.insert(0, '')
e_rvo2.bind("<Return>", (lambda event: update_changes(refresh=2)))
e_rvo2.pack(side='left', padx=20, pady=3, expand=1, fill='x')


#e_nset = Entry(tab_rv)
#e_nset.place(x=x0+xss, y=(sl+0)*lh, width=100)

e_sort = Entry(tab_rv)
e_sort.insert(0, 'BJD')
e_sort.place(x=x0+xss+sp*1, y=(12)*lh, width=100)

e_offset = Entry(tab_rv)
e_offset.insert(0, 400)
e_offset.bind("<Return>", (lambda event: update_changes()))
e_offset.place(x=x0+xss+sp*2, y=(12)*lh, width=100)

e_out = Entry(tab_rv)
e_out.insert(0, 'tmp.dat')
e_out.place(x=x0+xss, y=(13.5)*lh, width=200)


e_dir = Entry(tab_res)
e_dir.insert(0, 'res')
#e_dir.bind("<Return>", (lambda event: sel_res))
e_dir.place(x=x0+xss, y=2*lh, width=200)

e_nset_r = Entry(tab_res)
#e_nset_r.insert(0, '')
e_nset_r.place(x=x0+xss, y=3*lh, width=200)

e_oset_r = Entry(tab_res)
#e_oset_r.insert(0, '')
e_oset_r.place(x=x0+xss, y=4*lh, width=200)

###### COMBOBOXES ######

combo_avg = ttk.Combobox(tab_rv, values=['mean', 'wmean'])
combo_avg.set('wmean')
combo_avg.bind('<<ComboboxSelected>>', (lambda event: update_changes()))
combo_avg.place(x=x0+xss+sp*0, y=(12)*lh, width=100)


###### CHECKBOXES ######

ttk.Style().configure("TCheckbutton", background=bg_frame, bd=0, highlightthickness=0)

cb_ocen = IntVar()
ttk.Checkbutton(tab_rv, text="     ocen rvo 1", variable=cb_ocen).place(x=x0+sp, y=11*lh)

cb_cmpocen = IntVar()
ttk.Checkbutton(tab_rv, text="     ocen rvo 2", variable=cb_cmpocen).place(x=x0+sp*2, y=11*lh)

cb_cen = IntVar()
ttk.Checkbutton(tab_rv, text="     cen RVs to zero median", variable=cb_cen).place(x=x0, y=(11)*lh)

cb_cen.trace("w", update_changes)
cb_ocen.trace("w", update_changes)
cb_cmpocen.trace("w", update_changes)

###### LABELS ######

Label(master=win, text='Current command:').pack(side='top', anchor="w", padx=5)
e_run.pack(side='left', padx=5, pady=(0, 5))

Label(tab_rv, text='rvo file 1', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0, y=lh)

Label(tab_rv, text='rvo file 2', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0+420,y=lh)


Label(tab_rv, text='oset rvo 1:', background=bg_frame).place(x=x0, y=4*lh)
Label(tab_rv, text='oset rvo 2:', background=bg_frame).place(x=x0, y=7*lh)

#Label(tab_rv, text='Plot RV', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0, y=(sl-1.2)*lh)

#Label(tab_rv, text='nset:', background=bg_frame).place(x=x0, y=(sl+0)*lh)
Label(tab_rv, text='average:', background=bg_frame).place(x=x0+sp*0, y=(12)*lh)
Label(tab_rv, text='sort by:', background=bg_frame).place(x=x0+sp*1, y=(12)*lh)
Label(tab_rv, text='offset rvo:', background=bg_frame).place(x=x0+sp*2, y=(12)*lh)
Label(tab_rv, text='output:', background=bg_frame).place(x=x0, y=(13.5)*lh)

Label(tab_res, text='Plot residual', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0, y=lh)

Label(tab_res, text='directory:', background=bg_frame).place(x=x0, y=2*lh)
Label(tab_res, text='nset:', background=bg_frame).place(x=x0, y=3*lh)
Label(tab_res, text='oset:', background=bg_frame).place(x=x0, y=4*lh)


###### MAIN ######
if sys.argv[1:]:
    e_rvo1.delete(0, END)
    e_rvo1.insert(0, sys.argv[1])
refresh_oset('12')

win.mainloop()
