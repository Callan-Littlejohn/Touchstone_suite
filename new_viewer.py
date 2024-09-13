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
import darkdetect
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib
import configparser
import csv
import os
import scipy.optimize
from scipy.spatial.distance import cdist

class gui():
    def __init__(self):
        self.user_preferences={"bg":"","axes":"mz"}
        self.axes=self.user_preferences["axes"]
        self.figsize=(10,10)
        self.make_mainscreen()
    
    def make_mainscreen(self):
        maxrows=20
        maxcolumns=21
        self.mainscreen=tk.Tk()
        self.mainscreen.title("Touchstone")
        #self.mainscreen.configure(bg="DarkSlategrey")
        self.linedirection=tk.IntVar()
        self.mainscreen.iconbitmap("logo/touchstone_v2icon.ico")
        self.mainscreen.resizable(True,True)
        icon=tk.PhotoImage(file="logo/touchstone v2.png")
        self.mainscreen.iconphoto(False,icon )
        for i in range(maxrows):
            self.mainscreen.rowconfigure(i,weight=1)
        for i in range(maxcolumns):
            self.mainscreen.columnconfigure(i,weight=1)
        
        ## icons ##
        openpic=tk.PhotoImage(file="resources/pictures/newopenfile.png")
        openlabel=tk.Label(image=openpic)
        openbutton=tk.Button(self.mainscreen,image=openpic,command=self.open_file).grid(row=0,column=0)
        savepic=tk.PhotoImage(file="resources/pictures/newsavefile.png")
        savelabel=tk.Label(image=savepic)
        savebutton=tk.Button(self.mainscreen,image=savepic).grid(row=0,column=1)
        axespic=tk.PhotoImage(file="resources/pictures/axesicon.png")
        axeslabel=tk.Label(image=axespic)
        axesbutton=tk.Button(self.mainscreen,image=axespic).grid(row=0,column=2)
        peakpic=tk.PhotoImage(file="resources/pictures/pickpeakicon.png")
        peaklabel=tk.Label(image=peakpic)
        peakbutton=tk.Button(self.mainscreen,image=peakpic).grid(row=0,column=3)
        
        autocorpic=tk.PhotoImage(file="resources/pictures/autocor_icon.png")
        autocorlabel=tk.Label(image=autocorpic)
        autocorbutton=tk.Button(self.mainscreen,image=autocorpic,command=self.extract_ac).grid(row=0,column=15)
        
        horizpic=tk.PhotoImage(file="resources/pictures/horizontal_extract.png")
        horizlabel=tk.Label(image=horizpic)
        horizbutton=tk.Button(self.mainscreen,image=horizpic,command=self.show_horizontal).grid(row=0,column=16)
        
        vertpic=tk.PhotoImage(file="resources/pictures/vertical_extract.png")
        vertlabel=tk.Label(image=vertpic)
        vertbutton=tk.Button(self.mainscreen,image=vertpic).grid(row=0,column=17)
        self.lineentry=ttk.Entry(self.mainscreen)
        self.lineentry.grid(row=0,column=18,sticky="NSEW")
        self.lineentry.insert(0,"653.4")
        ## spectra tabs ##
        tabcontrol=ttk.Notebook(self.mainscreen)
        self.contourtab=ttk.Frame(tabcontrol)
        tabcontrol.add(self.contourtab,text="contour map") 
        self.spectrumtab=ttk.Frame(tabcontrol)
        tabcontrol.add(self.spectrumtab,text="spectrum")
        tabcontrol.grid(row=2,column=0,rowspan=15,columnspan=15,sticky="NESW")         
        x,y,z=np.asarray(range(100)),np.asarray(range(100)),[]
        for i in range(len(x)):
            z.append(y*i)
        self.contourfig=Figure(figsize=self.figsize)
        self.axcontour=self.contourfig.add_subplot(111)
        contour=self.axcontour.contourf(x,y,z)
        self.contourx,self.contoury=x,y
        self.contourcanvas=FigureCanvasTkAgg(self.contourfig,master=self.contourtab)
        self.contourcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.contourtoolbaer=NavigationToolbar2Tk(self.contourcanvas,self.contourtab)
        self.contourcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.contourcanvas.mpl_connect("button_press_event", self.on_contour_canvas_click)
        
        self.spectrumfig=Figure(figsize=self.figsize)
        self.axspectrum=self.spectrumfig.add_subplot(111)
        spectrum=self.axspectrum.plot(x,y)
        self.x,self.y=x,y
        self.spectrumcanvas=FigureCanvasTkAgg(self.spectrumfig,master=self.spectrumtab)
        self.spectrumcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.spectrumtoolbar=NavigationToolbar2Tk(self.spectrumcanvas,self.spectrumtab)
        self.spectrumcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.spectrumcanvas.mpl_connect("button_press_event", self.on_spectrum_canvas_click)
        
        ## mass lists ##
        tabcontrol2=ttk.Notebook(self.mainscreen)
        self.masslisttab=ttk.Frame(tabcontrol2)
        tabcontrol2.add(self.masslisttab,text="2D mass list")
        self.d1masslisttab=ttk.Frame(tabcontrol2)
        tabcontrol2.add(self.d1masslisttab,text="1D mass list")
        
        tabcontrol2.grid(row=2,column=15,columnspan=4,rowspan=15,sticky="NESW")
        
        self.masslistbox=ttk.Treeview(self.masslisttab)
        self.masslistbox["columns"]=("p_mass","f_mass","I")
        self.masslistbox.column("#0",width=0)
        self.masslistbox.column("p_mass",width=40)
        self.masslistbox.column("f_mass",width=40)
        self.masslistbox.column("I",width=40)
        self.masslistbox.heading("p_mass",text="P mass")
        self.masslistbox.heading("f_mass",text="F mass")
        self.masslistbox.heading("I",text="I")
        self.masslistbox.pack(fill=tk.BOTH,expand=1)
        
        self.d1masslistbox=ttk.Treeview(self.d1masslisttab)
        self.d1masslistbox["columns"]=("mass","I")
        self.d1masslistbox.column("#0",width=0)
        self.d1masslistbox.column("mass",width=40)
        self.d1masslistbox.column("I",width=40)
        self.d1masslistbox.heading("mass",text="mass")
        self.d1masslistbox.heading("I",text="I")
        self.d1masslistbox.pack(fill=tk.BOTH,expand=1)
        
        
        sv_ttk.set_theme("dark")
        self.mainscreen.mainloop()

    def open_file(self):
        self.thresholds=[5E6,5E7,7.5E7,1E8,1.5E8,5E8]
        dirname=tk.filedialog.askdirectory() # ask for file
        config=configparser.ConfigParser() #config parser tool
        config.read(str(dirname+"//parameters.p")) #opens config file
        filetype=config["spectrum_information"]["spectrum type"] # find filetype for future exploration, this allows one click opening
        self.data=np.abs(np.load(str(dirname+"//spectrum.npy")))
        print("dataimported")                
        with open(str(dirname+"//axesinfo.ax"),"r") as f:
            csv_reader=csv.reader(f, delimiter=",")
            axes=[]
            for row in csv_reader:
                axes.append(row)
            for i in range(len(axes)):
                for j in range(len(axes[i])):
                    if axes[i][j]=="":
                        axes[i][j]=0
                    axes[i][j]=float(axes[i][j])
        self.fragmentaxis=axes[0]
        self.precursoraxis=axes[1]
        for widget in self.contourtab.winfo_children():
            widget.destroy()
        self.contourfig=Figure(figsize=self.figsize)
        self.axcontour=self.contourfig.add_subplot(111)
        contourdata=np.where(self.data>self.thresholds[0])
        contourdata2 = (contourdata[0][self.data[contourdata] > 10], contourdata[1][self.data[contourdata] > 10])
        contoury=[self.precursoraxis[z] for z in contourdata2[0]]
        contourx=[self.fragmentaxis[z] for z in contourdata2[1]]
        contourz=self.data[contourdata2]
        cmin=min(contourz)
        cmax=max(contourz)
        contour = self.axcontour.scatter(contourx, contoury, c=contourz, cmap='plasma', s=np.divide(contourz,max(contourz)),norm=matplotlib.colors.LogNorm(vmin=cmin, vmax=cmax)) # Exchanging the contour plot for the scatter plot saves on time and processing power, this idea was taken from the PINK software produced by Anna Cordiner and all credit exists there
        self.contourx,self.contoury=self.fragmentaxis,self.precursoraxis
        self.axcontour.set_xlim([min(self.fragmentaxis[1:]),2000])
        self.axcontour.set_ylim([min(self.precursoraxis[1:]),2000])
        self.contourcanvas=FigureCanvasTkAgg(self.contourfig,master=self.contourtab)
        self.contourcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.contourtoolbaer=NavigationToolbar2Tk(self.contourcanvas,self.contourtab)
        self.contourcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.contourcanvas.mpl_connect("button_press_event", self.on_contour_canvas_click)
        self.contourcanvas.draw()
        self.clickpoint=0
        
    def save_file(self):
        dirname=tk.filedialog.asksaveasfilename() # ask for file
        dirname=dirname+".metal" #set filename with extension
        os.mkdir(dirname) #make the directory
        np.save(str(dirname+"/spectrum.npy"), self.data)
        with open(str(dirname+"/axesinfo.ax"),"w") as f: #writing axes
            f.write("0,")
            f.write(",".join([str(self.fragmentaxis[j]) for j in range(len(self.fragmentaxis))]))
            f.write(",\n")
            f.write("0,")
            f.write(",".join([str(self.precursoraxis[j]) for j in range(len(self.precursoraxis))]))
            f.close()
        configset=configparser.ConfigParser() 
        with open(str(dirname+"/parameters.p"),"w") as f: #setting up parameters
            configset["spectrum_information"]={"spectrum type":"2d","spectrum size": len(self.data[0]),"ysize":len(self.precursoraxis)}
            configset.write(f) #write params ths is a good sanity check
        print("saved")
    
    def on_spectrum_canvas_click(self, event):
        if event.inaxes:
            # Get the click coordinates
            print("yay")
            canvas_x, canvas_y = event.x, event.y
            # Convert canvas coordinates to graph data coordinates
            graph_x, graph_y = event.xdata, event.ydata
            # Get the nearest data point from the graph
            #closest_index = np.argmin(np.abs(self.x - graph_x))
            closest_index=np.argmin(np.abs(np.transpose(self.d1peaks)[0]-graph_x))
            closest_peak = self.d1peaks[closest_index]
            closest_data_x=closest_peak[0]
            closest_data_y=1.1*closest_peak[1]
            #self.last_spectrum_click=closest_index
            if hasattr(self,"mark"):
                if self.mark==0:
                    self.mark=self.axspectrum.scatter([closest_data_x],[closest_data_y],s=200,marker="v")
                self.mark.set_offsets(np.c_[closest_data_x,closest_data_y])
            else:
                self.mark=self.axspectrum.scatter([closest_data_x],[closest_data_y],s=200,marker="v")
            self.spectrumcanvas.draw()
            self.lineentry.delete(0,tk.END)
            self.lineentry.insert(0,str(closest_data_x))
            
    def on_contour_canvas_click(self, event):
        #print("yay")
        if event.inaxes:
            # Get the click coordinates
            canvas_x, canvas_y = event.x, event.y
            # Convert canvas coordinates to graph data coordinates
            graph_x, graph_y = event.xdata, event.ydata
            # Get the nearest data point from the graph
            closest_indexx = np.argmin(np.abs(self.contourx - graph_x))
            closest_indexy = np.argmin(np.abs(self.contoury - graph_y))
            closest_data_x = self.contourx[closest_indexx]
            closest_data_y = self.contoury[closest_indexy]
            self.last_contour_click_y,self.last_contour_click_x=closest_indexy,closest_indexx
            if hasattr(self,"clickpoint"):
                if self.clickpoint==0:
                    self.clickpoint=self.axcontour.scatter([closest_data_x],[closest_data_y],s=200,marker="x")
                    self.hline=self.axcontour.axhline(y=closest_data_y,c="r")
                    self.vline=self.axcontour.axvline(x=closest_data_x,c="r")
                self.clickpoint.set_offsets(np.c_[closest_data_x,closest_data_y])
                self.hline.set_ydata([closest_data_y,closest_data_y])
                self.vline.set_xdata([closest_data_x,closest_data_x])
            else:
                self.clickpoint=self.axcontour.scatter([closest_data_x],[closest_data_y],s=200,marker="x")
                self.hline=self.axcontour.axhline(y=closest_data_y,c="r")
                self.vline=self.axcontour.axvline(x=closest_data_x,c="r")
            self.contourcanvas.draw()
            self.lineentry.delete(0,tk.END)
            self.lineentry.insert(0,str(closest_data_y))
    def extract_ac(self):
        acxind=find_closest_indices(self.fragmentaxis,self.precursoraxis)
        ac=[]
        count=0
        for i in range(len(self.fragmentaxis)):
            ac.append(self.data[count,i])
            if i>acxind[count]:
                count=count+1
                if count==len(self.precursoraxis):
                    count=len(self.precursoraxis)-1
        x=self.fragmentaxis
        self.update_spectrum(x,ac)
    
    def show_horizontal(self):
        mzextract=float(self.lineentry.get())
        if mzextract==self.contoury[self.last_contour_click_y]:
            indextoshow=self.last_contour_click_y
        else:
            indextoshow=(np.abs(np.asarray(self.precursoraxis)-float(self.lineentry.get()))).argmin()-1
        x=self.fragmentaxis
        y=self.data[indextoshow]
        self.update_spectrum(x,y)
        
    def update_spectrum(self,x,y):
        for widget in self.spectrumtab.winfo_children():
            widget.destroy()
        peaksd=(peak_picker_guassian1d_wholespec([x[10:-10],y[10:-10]]))
        #print(peaks)
        self.spectrumfig=Figure(figsize=self.figsize)
        self.axspectrum=self.spectrumfig.add_subplot(111)
        spectrum=self.axspectrum.plot(x,np.abs(y),linewidth=0.5,color="black")
        if hasattr(self,"peaks")==False:
           self.peaksd1=peaksd
        else:
            self.peaksd1=peaksd
        self.d1peaks=self.peaksd1
        print(self.d1peaks)
        for i in self.d1masslistbox.get_children():
            self.d1masslistbox.delete(i)
        for i in self.peaksd1:
            try:
                self.d1masslistbox.insert("","end",str(i[1]),values=[i[0],i[1]])
            except:
                o=2
        self.axspectrum.set_xlim([min(x),2000])
        self.spectrumcanvas=FigureCanvasTkAgg(self.spectrumfig,master=self.spectrumtab)
        self.spectrumcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.spectrumtoolbar=NavigationToolbar2Tk(self.spectrumcanvas,self.spectrumtab)
        self.spectrumcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.spectrumcanvas.mpl_connect("button_press_event", self.on_spectrum_canvas_click)

def peak_picker_guassian1d_wholespec(spectrum):
    peaks=[]
    windowedge=3
    for i in range(30,len(spectrum[1])-31):
        if spectrum[0][i]>200 and spectrum[0][i]<3000:
            if spectrum[1][i]==max(spectrum[1][i-10:i+10]):
                if spectrum[1][i]>(np.mean(spectrum[1][i+3:i+10])+(2.5*(np.std(spectrum[1][i+3:i+10]))))*5:#(((np.sum(spectrum[1][i-15:i-2])-spectrum[1][i])/29)+(np.std(spectrum[1][i-15:i-1])*2))*5:
                    try:
                        xpeak=fitgaussian(spectrum[0][i-10:i+10], spectrum[1][i-10:i+10],po=[spectrum[1][i],spectrum[0][i],0.02])[1]
                        peaks.append([float(abs(xpeak)),spectrum[1][i]])
                    except:
                        print("couldnt work it out for:",spectrum[0][i])
    return peaks   
def fitgaussian(x,y,po=None):
    params,covar=scipy.optimize.curve_fit(gaussian,x,y,po)
    return params
def gaussian(x, amp, cent,width):
    return amp *np.exp(-((x-cent)**2)/(2*width**2))


def find_closest_indices(A, B):
    # Reshape arrays to make broadcasting work
    A=np.array(A)
    B=np.array(B)
    A_reshaped = A.reshape(-1, 1)
    B_reshaped = B.reshape(-1, 1)   
    differences = cdist(A_reshaped, B_reshaped) # Calculate absolute differences between A and B using cdist
    closest_indices = np.argmin(differences, axis=0)# Find indices with minimum differences for each point in B

    return closest_indices
gui()