#! /usr/bin/env python3
# -*- coding: iso-8859-1 -*-

# GUI to run vpr.py, showing most important options

import os
from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText
from hbox import Help_Box

import numpy as np
import sys
import vpr
vpr.pause = print   # mainloop of the gui will pause

viperdir = os.path.dirname(os.path.realpath(__file__)) + os.sep

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
        args = options[opt] + ' -save ' + e_out.get()
    elif '-plot par' in str(args):
        opt = '-par'
    else:
        opt = options.index(str(args))

    if not res and '-plot par' not in str(args):
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

    elif res:
        str_arg = " -res "+str(e_dir.get())

        if cb_resn != []:
            str_arg += " -nset ["
            for i,n in enumerate (cb_resn):
                if n.get():
                    print(i,n.get(),res_num)
                    str_arg += str(res_num[i])+","
            str_arg += "]"
        if cb_reso != []:
            str_arg += " -oset ["
            for i,o in enumerate (cb_reso):
                if o.get():
                    str_arg += str(res_orders[i])+","
            str_arg += "]"
        if e_ressep.get():
            str_arg += " -ressep "+e_ressep.get()

    #    if e_oset_r.get():
     #       str_arg += " -oset "+e_oset_r.get()
      #  if e_nset_r.get():
       #     str_arg += " -nset "+e_nset_r.get()
    else:
        str_arg = e_parf.get().split('.par')[0] 

    if args:
        str_arg += " " + args

    if '-plot par' in str(args): 
        str_arg += '-parcolx '+str(cb_parx.get())
        str_arg += ' -parcoly '+str(cb_pary.get())

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

def bt_par(e_file):
    file = askopenfilename()
    if file:
        e_file.delete(0, END)
        e_file.insert(0, file)

def bt_swap():
    f1 = filename1.get()
    f2 = filename2.get()
    filename1.set(f2)
    filename2.set(f1)
    refresh_oset('12')
    update_changes()

def text_from_file(text):
     with open(viperdir+"vpr.py") as search:
         for line in search:
             line = line.rstrip()  
             if text in line:
                  infotext = (line.split('help')[1]).split("'")[1]
     return infotext


cbo1, cbo2 = [], []			# oset checkboxes
cb_orders1, cb_orders2 = [], []		# Var() of  oset checkboxes 
o_rvo1, o_rvo2 = [], []			# orders in RV files
cb_reso, cb_resn = [], []
res1, res2 = [], []
res_num, res_orders = [], []

def refresh_oset(num):
    # refresh oset checkboxes (when using new RV file)
    if '1' in num:
        if e_rvo1.get():
            global cb_orders1, cbo1, o_rvo1 
            o_rvo1 = get_orders(e_rvo1.get())
            cb_orders1, cbo1 = create_cb(o_rvo1, cbo1, cb_orders1, lfr_oset1)
    if '2' in num:
        if e_rvo2.get():
            global cb_orders2, cbo2, o_rvo2 
            o_rvo2 = get_orders(e_rvo2.get())
            cb_orders2, cbo2 = create_cb(o_rvo2, cbo2, cb_orders2, lfr_oset2)

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

def create_cb(orders_all, cbo, cb_orders, frm):
    # creates checkboxes for oset
    if cbo != []:
        # delete old checkboxes
        for i in cbo:
            i.destroy()
        cb_orders = []
        cbo = []    

    for i, o in enumerate(orders_all):
        cb_orders.append(IntVar())
        yi, xi =  divmod(i, 5)
        c = ttk.Checkbutton(frm, text="  "+str(o),  variable=cb_orders[i])
        c.grid(row=yi, column=xi, sticky="nw", padx=15, pady=2)
        cbo.append(c)
        cb_orders[i].set(1)
        cb_orders[i].trace("w", update_changes)

    return cb_orders, cbo


def get_parameters():
    parfile = e_parf.get()#.split('.')[0]

    par = np.genfromtxt(parfile, dtype=None, names=True,
                            deletechars='',   # to keep the dash for chunks
                            encoding=None).view(np.recarray)
    colnames = par.dtype.names[:-1]

    s = ttk.Style()         
    s.configure('par.TRadiobutton', background=bg_frame)

    global cb_parx, cb_pary
    cb_parx, cb_pary = StringVar(), StringVar()

    # parameters x axis
    l = LabelFrame(frm_par1, text='', bg=bg_frame, bd=2)
    l.grid(row=1, column=0, sticky="news", padx=5, pady=5, ipady=5, columnspan=1)

    for i, c in enumerate(colnames[:3]):
        cc = ttk.Radiobutton(l, text="  "+str(c),  variable=cb_parx, value=c, style = 'par.TRadiobutton')
        cc.grid(row=i, column=0, sticky="new", padx=15, pady=2)
    cb_parx.set('n')

    # parameters y axis
    lfr_par = []
    par_groups = ['rv', 'norm', 'wave', 'ip', 'atm', 'bkg']
    for i, cn in enumerate(par_groups):
        yi, xi =  divmod(i, 6)
        l = LabelFrame(frm_par1, text=cn, bg=bg_frame, bd=2)
        l.grid(row=3+yi, column=xi, sticky="news", padx=5, pady=5, ipady=5, columnspan=1)
        frm_par1.grid_columnconfigure(i, weight=1)
        lfr_par.append(l)

    for i, c in enumerate(colnames[4:][::2]):
        pos = par_groups.index(re.sub(r'[0-9]', '', c))
        cc = ttk.Radiobutton(lfr_par[pos], text="  "+str(c),  variable=cb_pary, value=c, style = 'par.TRadiobutton')
        cc.grid(row=i+5, column=0, sticky="new", padx=15, pady=2)
    
    cb_parx.trace("w", lambda *args: call_vpr('-plot par '))
    cb_pary.trace("w", lambda *args: call_vpr('-plot par '))

    return cb_pary, cb_parx

def sel_res(cb_reso, cb_resn, res1, res2):
    list_dat = np.array(sorted(os.listdir(e_dir.get())))
    num, orders = [], []
    for dat in list_dat:
        no = (dat.split(".")[0]).split("_")
        num.append(int(no[0]))
        orders.append(int(no[1]))

    global res_num, res_orders
    res_num = np.unique(num)
    res_orders = np.unique(orders)

    for i, o in enumerate(res_orders):
        cb_reso.append(IntVar())
        yi, xi =  divmod(i, 5)
        c = ttk.Checkbutton(lreso, text="  "+str(o),  variable=cb_reso[i])
        c.grid(row=yi+1, column=xi, sticky="nw", padx=15, pady=2)
        res1.append(c)
        cb_reso[i].set(1)
        cb_reso[i].trace("w", lambda *args: call_vpr(res=True))
        lreso.grid_columnconfigure(xi, weight=1)

    for i, n in enumerate(res_num):
        cb_resn.append(IntVar())
        yi, xi =  divmod(i, 5)
        c = ttk.Checkbutton(lresn, text="  "+str(n),  variable=cb_resn[i])
        c.grid(row=yi+1, column=xi, sticky="nw", padx=15, pady=2)
        res2.append(c)
        if i == 0: cb_resn[0].set(1)
        cb_resn[i].trace("w", lambda *args: call_vpr(res=True))
        lresn.grid_columnconfigure(xi, weight=1)

    return res_num, res_orders


def bt_exit():
    exit()


###### SETTING ######

font_type = 'Arial'
font_size = 12
bg_frame = '#f1f1f1'    # bg color big frame
bg_color = '#e6e1e1'    # bg color small frames

win_width = 840     # width of GUI window
win_high = 820      # height of GUI window

win = Tk()
win.title('Gui VPR')
win.geometry("{}x{}".format(win_width, win_high))
win.configure(background=ttk.Style().lookup('TFrame', 'background'))   # same bg-color for win and Frame

style = ttk.Style(win)
style.configure('lefttab.TNotebook', tabposition='nw')
style.configure('TNotebook.Tab', font=(font_type, font_size+2, 'bold'))

tabControl = ttk.Notebook(win, style='lefttab.TNotebook')

tab_rv = ttk.Frame(tabControl)
tab_par = ttk.Frame(tabControl)
tab_res = ttk.Frame(tabControl)

tabControl.add(tab_rv, text ='RV plots')
tabControl.add(tab_par, text ='Parameters')
tabControl.add(tab_res, text ='Residuals')
tabControl.pack(expand=1, fill="both", side='top')


###### POSITIONS ######

fr_high = win_high-170
xy0 = 10

###### FRAMES #######

# options: solid, sunken, ridge, groove, raised
# groove, ridge need bd=2 or more

def new(event):
    win_width = win.winfo_width()
    win_height = event.height
    frm_rv.config(width=win_width-20)
    frm_par.config(width=win_width-20)
    frm_res.config(width=win_width-20)

win.bind("<Configure>", new)
win.columnconfigure(0, weight=1)

# Frame of RV Tab
frm_rv = Frame(tab_rv, height=fr_high, width=win_width-20, bg=bg_frame, bd=2, relief='groove')
frm_rv.grid(row=0, column = 0, sticky="news",padx=10,pady=10)
frm_rv.grid_propagate(False)
frm_rv.columnconfigure(0, weight=1)
frm_rv.columnconfigure(1, weight=1)

# Frame reading in rvo files
frm_rv1 = Frame(frm_rv, bg=bg_frame, bd=2, relief='groove')
#frm_rv1.pack(padx=10, pady=(10,10), ipady = (3), fill='both', side='top')
frm_rv1.grid(row=0, column=0, sticky="news", padx=10, pady=10, ipady=5, columnspan=2)
frm_rv1.grid_columnconfigure(1, weight=1)

# Sub Frames for RV parameters
lfr_oset1 = LabelFrame(frm_rv, text="oset rvo 1", bg=bg_frame, bd=2)
lfr_oset1.grid(row=1, column=0, sticky="news", padx=(10,0), pady=5, ipady=5, columnspan=1)

lfr_oset2 = LabelFrame(frm_rv, text="oset rvo 2", bg=bg_frame, bd=2)
lfr_oset2.grid(row=1, column=1, sticky="news", padx=10, pady=5, ipady=5, columnspan=1)

lfr_cen = LabelFrame(frm_rv, text="center RV values", bg=bg_frame, bd=2)
lfr_cen.grid(row=2, column=0, sticky="news", padx=(10,0), pady=5, ipady=5, columnspan=1)

lfr_plot = LabelFrame(frm_rv, text="Plotting Options", bg=bg_frame, bd=2)
lfr_plot.grid(row=2, column=1, sticky="news", padx=10, pady=5, ipady=5, columnspan=1)

lfr_other = LabelFrame(frm_rv, text="Other", bg=bg_frame, bd=2)
lfr_other.grid(row=3, column=0, sticky="news", padx=(10,10), pady=5, ipady=5, columnspan=2)

#lfr_other.grid_columnconfigure(1, weight=1)
#lfr_other.grid_columnconfigure(2, weight=1)

for cc in range(0,5,1):
    lfr_oset1.grid_columnconfigure(cc, weight=1)
    lfr_oset2.grid_columnconfigure(cc, weight=1)


# Frame for Parameter plots
frm_par = Frame(tab_par, height=100, width=win_width-20, bg=bg_frame, bd=2, relief='groove')
frm_par.grid(row=0, column = 0, sticky="news",padx=10,pady=(10,0))
frm_par.grid_propagate(False)
frm_par.grid_columnconfigure(1, weight=1)

frm_par1 = Frame(tab_par, height=fr_high-100, width=win_width-20, bg=bg_frame, bd=2, relief='groove')
frm_par1.grid(row=1, column = 0, sticky="news",padx=10,pady=(0,10))
frm_par1.grid_propagate(False)

# Frame for Residual Tab
frm_res = Frame(tab_res, height=fr_high, width=win_width-20, bg=bg_frame, bd=2, relief='groove')
frm_res.grid(row=0, column = 0, sticky="news",padx=10,pady=10)
frm_res.grid_propagate(False)
frm_res.grid_columnconfigure(0, weight=1)
frm_res.grid_columnconfigure(1, weight=1)

lres = LabelFrame(frm_res, text='Input', bg=bg_frame, bd=2)
lres.grid(row=1, column=0, sticky="news", padx=10, pady=5, ipady=5, columnspan=2)
#lres.grid_propagate(False)
lres.grid_columnconfigure(1, weight=1)

lreso = LabelFrame(frm_res, text='oset', bg=bg_frame, bd=2)
lreso.grid(row=7, column=1, sticky="news", padx=(5,10), pady=5, ipady=5, columnspan=1)

lresn = LabelFrame(frm_res, text='nset', bg=bg_frame, bd=2)
lresn.grid(row=7, column=0, sticky="news", padx=(10,5), pady=5, ipady=5, columnspan=1)


###### BUTTONS ######

ttk.Style().configure("TButton", padding=2,   background="#fdfdfd", font=(font_type,font_size,'bold'), borderwidth =2)
ttk.Style().configure("B2.TButton", padding=1,   background="#fdfdfd", font=(font_type,font_size-1,''), borderwidth =1)

b_exit = ttk.Button(master=win, text='EXIT', command = bt_exit)
b_exit.pack(side='right', anchor="s", padx=20, pady = 20)

b_rvo1 = Button(frm_rv1, text='Search data file', command = bt_rvo1, background="#fdfdfd", width=15)
b_rvo1.grid(row=0, column=2, sticky="ne", padx=xy0, pady=(xy0,3))

b_rvo2 = Button(frm_rv1,text='Search data file', command = bt_rvo2, background="#fdfdfd", width=15)
b_rvo2.grid(row=1, column=2, sticky="ne", padx=xy0, pady=(3))

b_swap = Button(frm_rv1, text='swap', command = bt_swap, background="#fdfdfd", width=15)
b_swap.grid(row=2, column=1, sticky="ne", padx=xy0, pady = 3, rowspan=2)

b_save = ttk.Button(lfr_other, text='Save', command = lambda: call_vpr('-save'))
b_save.grid(row=1, column=2, sticky="nw", padx=xy0)

b_cmp = ttk.Button(frm_rv1, text='Compare', command = lambda: call_vpr(cmp=True), width=15)
b_cmp.grid(row=2, column=2, sticky="ne", padx=xy0, pady = 3)

conf = {'side': 'right', 'anchor': 's', 'padx': (0,20), 'pady': 20}

b_nrvo = ttk.Button(frm_rv, text='Plot no-rv', command = lambda: call_vpr('-plot nrvo'))
b_nrvo.pack(**conf)

b_rvo = ttk.Button(frm_rv, text='Plot o-rv', command = lambda: call_vpr('-plot rvo'))
b_rvo.pack(**conf)

b_rvbjd = ttk.Button(frm_rv, text='Plot BJD-RV', command = call_vpr)
b_rvbjd.pack(**conf)

b_fpar = Button(frm_par,text='Search data file', command = lambda: bt_par(e_parf), background="#fdfdfd", width=15)
b_fpar.grid(row=1, column=2, sticky="ne", padx=10, pady=0)

#b_cmp.pack(**{**conf, 'side': 'right', 'padx': 20})

b_oset1_a = ttk.Button(lfr_oset1, text='select all', style = 'B2.TButton', command = lambda: set_oset(cb_orders1, 1))
b_oset1_a.grid(row=8, column=0, sticky="sw", padx=xy0, pady=5, columnspan=2)
b_oset1_n = ttk.Button(lfr_oset1, text='select none', style = 'B2.TButton', command = lambda: set_oset(cb_orders1, 0))
b_oset1_n.grid(row=8, column=2, sticky="sw", padx=xy0, pady=5, columnspan=2)

b_oset2_a = ttk.Button(lfr_oset2, text='select all', style = 'B2.TButton', command = lambda: set_oset(cb_orders2, 1))
b_oset2_a.grid(row=8, column=0, sticky="sw", padx=xy0, pady=5, columnspan=2)
b_oset2_n = ttk.Button(lfr_oset2, text='select none', style = 'B2.TButton', command = lambda: set_oset(cb_orders2, 0))
b_oset2_n.grid(row=8, column=2, sticky="sw", padx=xy0, pady=5, columnspan=2)

b_res = ttk.Button(frm_res, text='Plot res', command = lambda: call_vpr(res=True))
#b_res.grid(row=12, column=1, sticky="se", padx=xy0, pady=(xy0,xy0))
b_res.pack(**conf)

b_reso_a = ttk.Button(lreso, text='select all', style = 'B2.TButton', command = lambda: set_oset(cb_reso, 1))
b_reso_a.grid(row=0, column=0, sticky="ne", padx=20, pady=5, columnspan=6)
b_reso_n = ttk.Button(lreso, text='select none', style = 'B2.TButton', command = lambda: set_oset(cb_reso, 0))
b_reso_n.grid(row=0, column=0, sticky="ne", padx=110, pady=5, columnspan=6)

b_resn_a = ttk.Button(lresn, text='select all', style = 'B2.TButton', command = lambda: set_oset(cb_resn, 1))
b_resn_a.grid(row=0, column=0, sticky="ne", padx=20, pady=5, columnspan=6)
b_resn_n = ttk.Button(lresn, text='select none', style = 'B2.TButton', command = lambda: set_oset(cb_resn, 0))
b_resn_n.grid(row=0, column=0, sticky="ne", padx=110, pady=5, columnspan=6)


###### ENTRIES ######

filename1 = StringVar()
filename2 = StringVar()
parfile = StringVar()

Label(master=win, text='Current command:').pack(side='top', anchor="w", padx=xy0)
e_run = ScrolledText(win, background='#f0f0f0')
e_run.pack(side='left', padx=xy0, pady=(0, xy0))

e_rvo1 = Entry(frm_rv1, textvariable=filename1, width=100)
e_rvo1.insert(0, 'tmp.rvo.dat')
e_rvo1.bind("<Return>", (lambda event: update_changes(refresh=1)))
e_rvo1.grid(row=0, column=1, sticky="sew", padx=xy0, pady=(5))

e_rvo2 = Entry(frm_rv1, textvariable=filename2, width=100)
e_rvo2.insert(0, '')
e_rvo2.bind("<Return>", (lambda event: update_changes(refresh=2)))
e_rvo2.grid(row=1, column=1, sticky="new", padx=xy0, pady=(5))

e_sort = Entry(lfr_plot, width=10)
e_sort.insert(0, 'BJD')
e_sort.grid(row=0, column=1, sticky="nw", padx=xy0, pady=(5,0))

e_offset = Entry(lfr_plot, width=10)
e_offset.insert(0, 400)
e_offset.bind("<Return>", (lambda event: update_changes()))
e_offset.grid(row=1, column=1, sticky="nw", padx=(xy0), pady=(5,0))

e_out = Entry(lfr_other, width=35)
e_out.insert(0, 'tmp.dat')
e_out.grid(row=1, column=1, sticky="nw", padx=xy0, pady=(5,0), columnspan=1)


e_parf = Entry(frm_par, textvariable=parfile, width=120)
e_parf.insert(0, 'tmp.par.dat')
e_parf.bind("<Return>", (lambda event: get_parameters()))
e_parf.grid(row=1, column=1, sticky="new", padx=xy0, pady=(5))


e_dir = Entry(lres, width = 60)
e_dir.insert(0, 'res')
e_dir.bind("<Return>", (lambda event: sel_res(cb_reso, cb_resn, res1, res2)))
e_dir.grid(row=1, column=1, sticky="new", padx=xy0, pady=3)

e_ressep = Entry(lres, width = 20)
e_ressep.insert(0, 5)
e_ressep.bind("<Return>", (lambda event:  call_vpr(res=True)))
e_ressep.grid(row=2, column=1, sticky="nw", padx=xy0, pady=3)
#, lambda *args: call_vpr(res=True)
#e_nset_r = Entry(frm_res)
#e_nset_r.insert(0, '')
#e_nset_r.grid(row=2, column=1, sticky="nw", padx=xy0, pady=xy0)

#e_oset_r = Entry(frm_res)
#e_oset_r.insert(0, '')
#e_oset_r.grid(row=3, column=1, sticky="nw", padx=xy0, pady=xy0)


###### COMBOBOXES ######

combo_avg = ttk.Combobox(lfr_other, values=['mean', 'wmean'], width=8)
combo_avg.set('wmean')
combo_avg.bind('<<ComboboxSelected>>', (lambda event: update_changes()))
combo_avg.grid(row=0, column=1, sticky="nw", padx=xy0, pady=(5,0))


###### CHECKBOXES ######

ttk.Style().configure("TCheckbutton", background=bg_frame, bd=0, highlightthickness=0)

cb_ocen = IntVar()
l_ocen1 = ttk.Checkbutton(lfr_cen, text="     ocen rvo 1", variable=cb_ocen)
l_ocen1.grid(row=1, column=0, sticky="nw", padx=xy0, pady=(5,0))
Help_Box(widget = l_ocen1, text = text_from_file("'-ocen'"))

cb_cmpocen = IntVar()
l_ocen2 = ttk.Checkbutton(lfr_cen, text="     ocen rvo 2", variable=cb_cmpocen)
l_ocen2.grid(row=2, column=0, sticky="nw", padx=xy0, pady=(5,0))
Help_Box(widget = l_ocen2, text = text_from_file("'-cmpocen'"))

cb_cen = IntVar()
l_cen = ttk.Checkbutton(lfr_cen, text="     cen RVs to zero median", variable=cb_cen)
l_cen.grid(row=0, column=0, sticky="nw", padx=xy0, pady=(5,0))
Help_Box(widget = l_cen, text = text_from_file("'-cen'"))

cb_cen.trace("w", update_changes)
cb_ocen.trace("w", update_changes)
cb_cmpocen.trace("w", update_changes)


###### LABELS ######

Label(frm_rv1, text='rvo file 1', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=0, column=0, sticky="sw", padx=xy0, pady=5)

Label(frm_rv1, text='rvo file 2', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=1, column=0, sticky="nw", padx=xy0, pady=5)

l_sort = Label(lfr_plot, text='sort by:', background=bg_frame)
l_sort.grid(row=0, column=0, sticky="nw", padx=xy0, pady=(5,0))
Help_Box(widget = l_sort, text = text_from_file("'-sort'"))

l_off = Label(lfr_plot, text='offset rvo:', background=bg_frame)
l_off.grid(row=1, column=0, sticky="nw", padx=xy0, pady=(5,0))
Help_Box(widget = l_off, text = text_from_file("'-offset'"))

l_aver = Label(lfr_other, text='average:', background=bg_frame)
l_aver.grid(row=0, column=0, sticky="nw", padx=xy0, pady=(5,0))
Help_Box(widget = l_aver, text = text_from_file("'-avg'"))

l_out = Label(lfr_other, text='output:', background=bg_frame)
l_out.grid(row=1, column=0, sticky="nw", padx=xy0, pady=(5,0))
Help_Box(widget = l_out, text = text_from_file("'-save'"))

Label(frm_par, text='Plot parameters', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=0, column=0, sticky="nw", padx=xy0, pady=xy0, columnspan = 2)

Label(frm_par, text='parfile:', background=bg_frame).grid(row=1, column=0, sticky="nw", padx=xy0, pady=5)

Label(frm_par1, text='Parameters x-axis', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=0, column=0, sticky="nw", padx=xy0, pady=xy0, columnspan = 5)

Label(frm_par1, text='Parameters y-axis', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=2, column=0, sticky="nw", padx=xy0, pady=xy0, columnspan = 5)

Label(frm_res, text='Plot residual', font=(font_type, font_size, 'bold'), background=bg_frame).grid(row=0, column=0, sticky="nw", padx=xy0, pady=xy0)

Label(frm_res, text='For large number of nset, it is recommended to first chose one single order before using "select all" on nset.', background=bg_frame, fg='red4').grid(row=2, column=0, sticky="nw", padx=xy0, pady=xy0, columnspan=2)

Label(lres, text='directory:', background=bg_frame).grid(row=1, column=0, sticky="nw", padx=xy0, pady=3)

l_ressep = Label(lres, text='offset res:', background=bg_frame)
l_ressep.grid(row=2, column=0, sticky="nw", padx=xy0, pady=3)
Help_Box(widget = l_ressep, text = text_from_file("'-ressep'"))
#Label(frm_res, text='nset:', background=bg_frame).grid(row=2, column=0, sticky="nw", padx=xy0, pady=xy0)
#Label(frm_res, text='oset:', background=bg_frame).grid(row=3, column=0, sticky="nw", padx=xy0, pady=xy0)


###### MAIN ######
if sys.argv[1:]:
    e_rvo1.delete(0, END)
    e_rvo1.insert(0, sys.argv[1])

if os.path.isfile(e_rvo1.get()):
    refresh_oset('12')

if os.path.isfile(e_parf.get()):
    get_parameters()

if os.path.isdir(e_dir.get()):
    sel_res(cb_reso, cb_resn, res1, res2)

win.mainloop()
