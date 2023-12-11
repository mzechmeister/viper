#! /usr/bin/env python3

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

import tkinter as tk

#infotext
class Help_Box(object):
    def __init__(self, widget, text, show=1):
        self.widget = widget
        self.text = text

        def show_help(event):
            self.iwindow = tw = tk.Toplevel(self.widget)
            tw.wm_overrideredirect(1) 
            tw.wm_geometry("+{}+{}".format(self.widget.winfo_rootx()+30, self.widget.winfo_rooty()+20))
            label = tk.Label(tw, text = self.text, background = "#ffffe0", relief = 'solid', borderwidth = 1, wraplength=300, justify="left").pack(ipadx=10, ipady=2)

        def no_help(event):
            tw = self.iwindow
            tw.destroy()
            self.iwindow = None

        if show:
            widget.bind('<Enter>', show_help)
            widget.bind('<Leave>', no_help)
