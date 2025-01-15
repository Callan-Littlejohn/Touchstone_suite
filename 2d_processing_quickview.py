# -*- coding: utf-8 -*-
"""
----License----

Copyright Ó 2024  Callan Littlejohn,  Peter O'Connor, The University of Warwick and Verdel instruments

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or  any later version. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

For more details, see the GNU General Public License.

 

----Notes----

This software was produced by Callan Littlejohn as part of a collaboration between Verdel Instruments and the University of Warwick. It is intended to make the processing and analysis of 2DMS data simpler and faster.

"""

import tkinter as tk
from tkinter import ttk
import sv_ttk
from twoDMS_processing import d2spectrum
from tkinter import filedialog
import time

class d2msquickviewgui():
    def __init__(self):      
        self.mainscreen=tk.Tk()
        self.filetype=tk.StringVar()
        self.mainscreen.title("2DMS Processing Quickview")
        self.welcomelabel=ttk.Label(self.mainscreen, text="2DMS processing")
        self.ready=tk.Label(self.mainscreen,bg="green")
        self.ready.grid(row=1,column=0,rowspan=3,sticky="nsew")
        self.infilelabel=ttk.Label(self.mainscreen,text="please choose an infile")
        self.infilelabel.grid(row=1,column=1)
        self.outfilelabel=ttk.Label(self.mainscreen,text="please choose a destination")
        self.outfilelabel.grid(row=2,column=1)
        self.namelabel=ttk.Label(self.mainscreen,text="name").grid(row=3,column=1)
        
        self.infilebrowse=ttk.Button(self.mainscreen,text="browse",command=self.browseinfile).grid(row=1,column=2,sticky="NSEW")
        self.outfilebrowse=ttk.Button(self.mainscreen,text="browse",command=self.browseoutfile).grid(row=2,column=2,sticky="NSEW")
        self.nameentry=ttk.Entry(self.mainscreen)
        self.nameentry.grid(row=3,column=2)
        
        self.csvradio=ttk.Radiobutton(self.mainscreen,text="csv",variable=self.filetype,value="csv")
        self.npradio=ttk.Radiobutton(self.mainscreen,text="npy",variable=self.filetype,value="npy")
        self.csvradio.grid(row=4,column=2)
        self.npradio.grid(row=4,column=1)
        self.npradio.invoke()
        self.gobutton=ttk.Button(self.mainscreen,text="GO!",command=self.gotime).grid(row=1,column=3,rowspan=3,sticky="nsew")
        
        sv_ttk.set_theme("dark")
        self.mainscreen.mainloop()
        
    def browseinfile(self):
        self.infile=filedialog.askdirectory()
        self.infilelabel.configure(text=self.infile)
    
    def browseoutfile(self):
        self.outfile=filedialog.askdirectory()
        self.outfilelabel.configure(text=self.outfile)
    
    def gotime(self):
        t0=time.time()
        self.ready.config(bg="red")
        name=self.nameentry.get()
        finalout=self.outfile+"/"+name
        d2class=d2spectrum(self.infile)
        d2class.process2d()
        d2class.save2d(finalout,filetype=self.filetype.get())
        self.ready.config(bg="green")
        t1=time.time()
        print("done:",t1-t0)
    
if __name__=="__main__":
    r=d2msquickviewgui()