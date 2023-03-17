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
