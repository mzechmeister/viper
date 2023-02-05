#! /usr/bin/env python3
# -*- coding: iso-8859-1 -*-

# GUI to run vpr.py, showing most important options

import os
from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter import ttk

import sys
import vpr
vpr.pause = print   # mainloop of the gui will pause

###### FUNCIONS BUTTONS ######

def bt_start(args='-plot rv', cmp=False, res=False):
    print('---Starting vpr.py with choosen parameters---')

    if not res:
        str_arg = e_rvo1.get()

        if cb_cen.get():
            str_arg += " -cen "
        if cb_ocen.get():
            str_arg += " -ocen "
        if e_oset1.get():
            str_arg += " -oset "+e_oset1.get()
        # if e_nset.get():
        #     str_arg += " -nset "+str(e_nset.get())
        if e_sort.get():
            str_arg += " -sort "+str(e_sort.get())

        if cmp:
            if e_rvo2.get():
               str_arg += " -cmp "+e_rvo2.get()
            if e_oset2.get():
                str_arg += " -cmposet "+e_oset2.get()
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

    print("option string:", str_arg)
    vpr.run(str_arg.split())

    print('---Finished vpr.py---')

def bt_rvo1():
    file1 = askopenfilename(master=win)
    if file1:
        e_rvo1.delete(0, END)
        e_rvo1.insert(0, file1)

def bt_rvo2():
    file2 = askopenfilename(master=win)
    if file2:
        e_rvo2.delete(0, END)
        e_rvo2.insert(0, file2)

def bt_exit():
    exit()


###### SETTING ######

font_type = 'Arial'
font_size = 12
bg_frame = '#f1f1f1'    # bg color big frame
bg_color = '#e6e1e1'    # bg color small frames

win_width = 840     # width of GUI window
win_high = 405      # height of GUI window

win = Tk()
win.title('Gui VPR')
win.geometry("{}x{}".format(win_width, win_high))
win.configure(background=bg_color)


###### POSITIONS ######

x0 = 40     # positions start left
lh = 30     # line height
xs = 100
xss = 70
sp = 220
sps = 170
sl = 7

###### FRAMES #######

# options: solid, sunken, ridge, groove, raised
# groove, ridge need bd=2 or more
Frame(master=win, height=138, width=win_width-40, bg=bg_frame, bd=2, relief='groove').place(x=20, y=20)

Frame(master=win, height=185, width=390, bg=bg_frame, bd=2, relief='groove').place(x=20, y=50+112)
Frame(master=win, height=185, width=404, bg=bg_frame, bd=2, relief='groove').place(x=415, y=50+112)


###### BUTTONS ######

ttk.Style().configure("TButton", padding=2,   background="#fdfdfd", font=(font_type,font_size,'bold'), borderwidth =2)

b_exit = ttk.Button(master=win, text='EXIT', command = bt_exit)
b_exit.place(x=win_width-82,y=win_high-50, width=60)

b_rvo1 = Button(master=win, text='Search data file', command = bt_rvo1, background="#fdfdfd")
b_rvo1.place(x=win_width-580, y=lh, height=25)

b_rvo2 = Button(master=win,text='Search data file', command = bt_rvo2, background="#fdfdfd")
b_rvo2.place(x=win_width-179, y=lh, height=25)

call_vpr = bt_start
width_bt = 110

b_rvbjd = ttk.Button(master=win, text='Plot BJD-RV', command = call_vpr)
b_rvbjd.place(x=x0, y=(sl+3)*lh, width=width_bt)

b_rvo = ttk.Button(master=win, text='Plot o-rv', command = lambda: call_vpr('-plot rvo'))
b_rvo.place(x=x0+(width_bt+10)*1, y=(sl+3)*lh, width=width_bt)

b_cmp = ttk.Button(master=win, text='Compare', command = lambda: call_vpr(cmp=True))
b_cmp.place(x=x0+(width_bt+10)*2, y=(sl+3)*lh, width=width_bt)

b_res = ttk.Button(master=win, text='Plot res', command = lambda: call_vpr(res=True))
b_res.place(x=win_width-width_bt-x0, y=(sl+3)*lh, width=width_bt)

###### ENTRIES ######

e_rvo1 = Entry(master=win)
e_rvo1.insert(0, 'tmp.rvo.dat')
e_rvo1.place(x=x0, y=2*lh, width=350)

e_rvo2 = Entry(master=win)
e_rvo2.insert(0, '')
e_rvo2.place(x=x0+400, y=2*lh, width=350)


e_oset1 = Entry(master=win)
e_oset1.place(x=x0+xs, y=3*lh, width=150)

e_oset2 = Entry(master=win)
e_oset2.place(x=x0+xs+400, y=3*lh, width=150)

#e_nset = Entry(master=win)
#e_nset.place(x=x0+xss, y=(sl+0)*lh, width=100)

e_sort = Entry(master=win)
e_sort.insert(0, 'BJD')
e_sort.place(x=x0+xss, y=(sl+1)*lh, width=100)


e_dir = Entry(master=win)
e_dir.insert(0, 'res/')
e_dir.place(x=x0+xss+400, y=(sl+0)*lh, width=100)

e_nset_r = Entry(master=win)
#e_nset_r.insert(0, '')
e_nset_r.place(x=x0+xss+400, y=(sl+1)*lh, width=100)

e_oset_r = Entry(master=win)
#e_oset_r.insert(0, '')
e_oset_r.place(x=x0+xss+400, y=(sl+2)*lh, width=100)

###### COMBOBOXES ######


###### CHECKBOXES ######

ttk.Style().configure("TCheckbutton", background=bg_frame, bd=0, highlightthickness=0)

cb_ocen = IntVar()
ttk.Checkbutton(master=win, text="     center orders rvo 1", variable=cb_ocen).place(x=x0, y=4*lh)

cb_cmpocen = IntVar()
ttk.Checkbutton(master=win, text="     center orders rvo 2", variable=cb_cmpocen).place(x=x0+400, y=4*lh)

cb_cen = IntVar()
ttk.Checkbutton(master=win, text="     center RVs to zero median", variable=cb_cen).place(x=x0, y=(sl+0)*lh)

###### LABELS ######

Label(master=win, text='rvo file 1', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0, y=lh)

Label(master=win, text='rvo file 2', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0+400,y=lh)


Label(master=win, text='oset rvo 1:', background=bg_frame).place(x=x0, y=3*lh)
Label(master=win, text='oset rvo 2:', background=bg_frame).place(x=x0+400, y=3*lh)

Label(master=win, text='Plot RV', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0, y=(sl-1.2)*lh)

Label(master=win, text='Plot residual', font=(font_type, font_size, 'bold'), background=bg_frame).place(x=x0+sp+sps, y=(sl-1.2)*lh)

#Label(master=win, text='nset:', background=bg_frame).place(x=x0, y=(sl+0)*lh)
Label(master=win, text='sort by:', background=bg_frame).place(x=x0, y=(sl+1)*lh)

Label(master=win, text='directory:', background=bg_frame).place(x=x0+sp+sps, y=(sl+0)*lh)
Label(master=win, text='nset:', background=bg_frame).place(x=x0+sp+sps, y=(sl+1)*lh)
Label(master=win, text='oset:', background=bg_frame).place(x=x0+sp+sps, y=(sl+2)*lh)


###### MAIN ######
if sys.argv[1:]:
    e_rvo1.delete(0, END)
    e_rvo1.insert(0, sys.argv[1])

win.mainloop()
