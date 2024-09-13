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
import matplotlib.pyplot as plot
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
import configparser
import csv
import scipy
import scipy.optimize
import scipy.signal
from scipy import optimize
import os
import math
from scipy.spatial.distance import cdist
import matplotlib
class showinggui():
    def __init__(self):
        back="darkblue"
        fore="azure3"
        self.peaks=[]
        self.figsize=(10,10)
        ## initialisation##        
        self.mainscreen=tk.Tk()
        self.mainscreen.title("2DMS viewer")
        self.mainscreen.configure(bg=back)
        self.linedirection=tk.IntVar()
        ## top bar ##
        topbarmenu=tk.Menu(self.mainscreen,activebackground="Ivory4",bg="Ivory3")
        self.mainscreen.config(menu=topbarmenu)
        filebar=tk.Menu(topbarmenu)
        editbar=tk.Menu(topbarmenu)
        toolbar=tk.Menu(topbarmenu)
        
        filebar.add_command(label="open", command=self.open_file)
        filebar.add_command(label="save", command=self.save_file)
        
        toolbar.add_command(label="find peaks",command=self.findpeaks2dguass)
        toolbar.add_command(label="calibration",command=self.calibration_tab)
        
        topbarmenu.add_cascade(label="File",menu=filebar)
        topbarmenu.add_cascade(label="Edit",menu=editbar)
        topbarmenu.add_cascade(label="Tools",menu=toolbar)
        ## icons ##
        openpic=tk.PhotoImage(file="resources/pictures/openfile.png")
        openlabel=tk.Label(image=openpic)
        openbutton=tk.Button(self.mainscreen,image=openpic,command=self.open_file).grid(row=0,column=1)
        savepic=tk.PhotoImage(file="resources/pictures/savefile.png")
        savelabel=tk.Label(image=savepic)
        savebutton=tk.Button(self.mainscreen,image=savepic).grid(row=0,column=2)
        horizontalselector=tk.Radiobutton(self.mainscreen,text="horizontal",padx=20,variable=self.linedirection,value=1,fg=fore,bg=back,selectcolor="black")
        horizontalselector.grid(row=0,column=5)
        verticalselector=tk.Radiobutton(self.mainscreen,text="vertical",padx=20,variable=self.linedirection,value=2,fg=fore,bg=back,selectcolor="black")
        verticalselector.grid(row=0,column=6)
        acselector=tk.Radiobutton(self.mainscreen,text="autocorr",padx=20,variable=self.linedirection,value=3,fg=fore,bg=back,selectcolor="black")
        acselector.grid(row=0,column=7)
        acselector.select()
        kmdlabel=tk.Label(self.mainscreen,text="kmd ref mass",fg=fore,bg=back).grid(row=0,column=10)
        self.kmdrefentry=tk.Entry(self.mainscreen,bd=5)
        self.kmdrefentry.grid(row=0,column=11,sticky="NESW")
        self.kmdrefentry.insert(0,"12.0")
        linefinderlabel=tk.Label(self.mainscreen,text="m/z",fg=fore,bg=back).grid(row=0,column=12)
        self.linefinderentry=tk.Entry(self.mainscreen,bd=5)
        self.linefinderentry.grid(row=0,column=13,columnspan=2,sticky="NESW")
        self.linefinderentry.insert(0,"653.4")
        lineselectorbutton=tk.Button(self.mainscreen,text="select", command=self.show_line,bg=back,fg=fore).grid(row=0,column=15,sticky="NESW")
        ## spectra tabs ##
        tabcontrol=ttk.Notebook(self.mainscreen)
        self.contourtab=ttk.Frame(tabcontrol)
        tabcontrol.add(self.contourtab,text="contour map") 
        self.spectrumtab=ttk.Frame(tabcontrol)
        tabcontrol.add(self.spectrumtab,text="spectrum")
        self.kmdtab=ttk.Frame(tabcontrol)
        tabcontrol.add(self.kmdtab,text="kmd")
        tabcontrol.grid(row=2,column=1,rowspan=15,columnspan=11,sticky="NESW")         
        x,y,z=np.asarray(range(100)),np.asarray(range(100)),[]
        for i in range(len(x)):
            z.append(y*i)
        self.contourfig=Figure(figsize=self.figsize)
        self.axcontour=self.contourfig.add_subplot(111)
        contour=self.axcontour.contourf(x,y,z)
        self.contourcanvas=FigureCanvasTkAgg(self.contourfig,master=self.contourtab)
        self.contourcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.contourtoolbaer=NavigationToolbar2Tk(self.contourcanvas,self.contourtab)
        self.contourcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        
        self.spectrumfig=Figure(figsize=self.figsize)
        self.axspectrum=self.spectrumfig.add_subplot(111)
        spectrum=self.axspectrum.plot(x,y)
        self.spectrumcanvas=FigureCanvasTkAgg(self.spectrumfig,master=self.spectrumtab)
        self.spectrumcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.spectrumtoolbar=NavigationToolbar2Tk(self.spectrumcanvas,self.spectrumtab)
        self.spectrumcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        
        self.kmdfig=Figure(figsize=self.figsize)
        self.axspectrum=self.kmdfig.add_subplot(111)
        kmd=self.axspectrum.scatter(x,y)
        self.kmdcanvas=FigureCanvasTkAgg(self.kmdfig,master=self.kmdtab)
        self.kmdcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.kmdtoolbar=NavigationToolbar2Tk(self.kmdcanvas,self.kmdtab)
        self.kmdcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        ## mass lists ##
        tabcontrol2=ttk.Notebook(self.mainscreen)
        self.masslisttab=ttk.Frame(tabcontrol2)
        tabcontrol2.add(self.masslisttab,text="2D mass list")
        self.d1masslisttab=ttk.Frame(tabcontrol2)
        tabcontrol2.add(self.d1masslisttab,text="1D mass list")
        
        tabcontrol2.grid(row=2,column=12,columnspan=4,rowspan=15,sticky="NESW")
        
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
        
        
        self.mainscreen.mainloop()
        
        
    def open_file(self):
        self.thresholds=[5E7,5E7,7.5E7,1E8,1.5E8,5E8]
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
        print(len(self.data))
        print(len(self.data[1]))
        print(len(axes[1]))
        print(len(axes[0]))
        self.fragmentaxis=axes[0]
        self.precursoraxis=axes[1]
        print(np.max(self.data))
        print(np.min(self.data))
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
        #contour=self.axcontour.imshow()
        self.axcontour.set_xlim([min(self.fragmentaxis[1:]),2000])
        self.axcontour.set_ylim([min(self.precursoraxis[1:]),2000])
        self.contourcanvas=FigureCanvasTkAgg(self.contourfig,master=self.contourtab)
        self.contourcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.contourtoolbaer=NavigationToolbar2Tk(self.contourcanvas,self.contourtab)
        self.contourcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.contourcanvas.draw()
        
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
    
    def filtering_vertical_ridges_in_1D(self,peaksx):
        print("filtering")
        yi=(np.abs(np.asarray(self.precursoraxis)-float(self.linefinderentry.get()))).argmin()
        ym=self.precursoraxis[yi]
        print(ym)
        peaks3dinline=[]
        peakys=np.transpose(self.peaks)[0]
        peakys=np.abs(np.subtract(peakys,ym))
        locator=np.where(peakys<5)[0]
        print(locator)
        peaksinx=[]
        peaksfinal=[]
        vertridgepeaks=[]
        for i in locator:
            peaksinx.append(self.peaks[i][1])
            print(self.peaks[i][0],self.peaks[i][1])
        for i in range(len(peaksx[0])):
            if min(np.subtract(peaksinx,peaksx[1][i]))<5:
                peaksfinal.append([peaksx[0][i],peaksx[1][i]])
            else:
                vertridgepeaks.append([peaksx[0][i],peaksx[1][i]])
        return np.transpose(peaksfinal),np.transpose(vertridgepeaks)      
        
    def show_line(self):
        print(self.linedirection.get())
        if self.linedirection.get()==1:
            index=(np.abs(np.asarray(self.precursoraxis)-float(self.linefinderentry.get()))).argmin()-1
            x=self.fragmentaxis
            y=self.data[index]
        if self.linedirection.get()==2:
            index=(np.abs(np.asarray(self.fragmentaxis)-float(self.linefinderentry.get()))).argmin()
            x=self.precursoraxis
            y=self.data[:,index]
        if self.linedirection.get()==3:
            x=self.fragmentaxis
            y=self.extract_ac()
        for widget in self.spectrumtab.winfo_children():
            widget.destroy()
        for widget in self.kmdtab.winfo_children():
            widget.destroy()
        self.linex=x
        self.liney=y
        peaksd=(peak_picker_guassian1d_wholespec([x[10:-10],y[10:-10]]))
        #print(peaks)
        self.spectrumfig=Figure(figsize=self.figsize)
        self.axspectrum=self.spectrumfig.add_subplot(111)
        spectrum=self.axspectrum.plot(x,np.abs(y),linewidth=0.5,color="black")
        if len(self.peaks)>0:
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
        ykmd=findkmdarray(float(self.kmdrefentry.get()),np.transpose(self.peaksd1))
        self.kmdfig=Figure(figsize=self.figsize)
        self.kmdspectrum=self.kmdfig.add_subplot(111)
        kmd=self.kmdspectrum.scatter(np.transpose(self.peaksd1)[0],ykmd)
        self.kmdspectrum.set_xlim([min(x),2000])
        self.kmdcanvas=FigureCanvasTkAgg(self.kmdfig,master=self.kmdtab)
        self.kmdcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.kmdtoolbar=NavigationToolbar2Tk(self.kmdcanvas,self.kmdtab)
        self.kmdcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)

        self.spectrumcanvas.draw()     
        
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
        return ac
        
    def locatenoisebins(self,xaxes,yaxes,binsizemzx,binsizemzy):
        x=1E10
        y=1E10
        xbins,ybins=[],[]
        xtemp,ytemp=[],[]
        xpos,ypos=[],[]
        xbincount,ybincount=[],[]
        for i in range(len(xaxes)):
            if xaxes[i]<3000:
                if abs(xaxes[i]-x)>binsizemzx:
                    xbins.append(i)
                    x=xaxes[i]
                    xtemp.append(xaxes[i])
                    xpos.append(xtemp)
                    xtemp=[]
                xtemp.append(xaxes[i])
        for i in range(len(yaxes)):
            if abs(yaxes[i]-y)>binsizemzy:
                ybins.append(i)
                y=yaxes[i]
                ytemp.append(xaxes[i])
                ypos.append(xtemp)
                ytemp=[]
            ytemp.append(xaxes[i])
        return xbins,ybins,xpos,ypos          
                
    def findpeaks2dbasic(self):
        self.peaks=peak_picker_basic_3d_optimised([self.fragmentaxis,self.precursoraxis,self.data])
        for i in self.peaks:
            self.masslistbox.insert("","end",str(i[2]),values=[i[0],i[1],i[2]])
        
    def findpeaks2dguass(self):
        self.peaks=peak_picker_guassian([self.fragmentaxis[10:],self.precursoraxis[10:],self.data[10:,10:]])
        for i in self.peaks:
            self.masslistbox.insert("","end",str(i[2]),values=[i[0],i[1],i[2]])
        
    def findpeaks2d(self):
        snratio=10
        noiselevels=[]
        positions=[]
        pos2d=[]
        peaklocations=[]
        peaks=[]
        pospicked2d=[]
        xbins,ybins,xpos,ypos,=self.locatenoisebins(self.fragmentaxis,self.precursoraxis,0.3,5)
        for i in range(1,len(ybins)-2):
            temp=[]
            pos2dtemp=[]
            for j in range(1,len(xbins)-2):
                temp.append(np.average(self.data[ybins[i]:ybins[i+1],xbins[j]:xbins[j+1]]))
                pos2dtemp.append([ypos[i],xpos[j]])
            pos2d.append(pos2dtemp)
            noiselevels.append(temp)
        noiselevels=np.asanyarray(noiselevels)
        for i in range(5,len(noiselevels)-5):
            for j in range(5,len(noiselevels[i]-5)):
                if noiselevels[i,j]>10*(np.average(noiselevels[i-5:i+5,j-5:j+5])):
                    peaklocations.append([i,j])
                    pospicked2d.append(pos2d[i][j])
        for i in peaklocations:
            locatory=ybins[i[0]]
            locatorx=xbins[i[1]]
            area=self.data[ybins[i[0]]:ybins[i[0]+1],xbins[i[1]]:xbins[i[1]+1]]
            yfinder=[]
            for j in area:
                yfinder.append(np.mean(j))
            yrough=np.argmax(yfinder)
            xpeak=fitgaussian(self.fragmentaxis[xbins[i[1]]],area[yrough],[1E7,self.fragmentaxis[xbins[i[1]]],1])
            peaks.append(xpeak)
            print(xpeak)
        print(peaks)
    
    def calibration_tab(self):
        back="darkblue"
        fore="azure3"        
        self.removed=[]
        self.callibration_tab=tk.Tk()
        self.callibration_tab.title("2DMS calibration")
        self.callibration_tab.configure(bg=back)
        self.callibration_lists=[f for f in os.listdir("reflists") if os.path.isfile(os.path.join("reflists",f))]
        self.callibration_combobox=ttk.Combobox(self.callibration_tab,values=self.callibration_lists)
        self.callibration_combobox.grid(row=0,column=0,columnspan=3)
        callibration_loadlist_button=tk.Button(self.callibration_tab,text="select list",command=self.select_cal_list).grid(row=0,column=3)
        self.calibration_box=ttk.Treeview(self.callibration_tab)
        self.calibration_box["columns"]=("Label","Ref Mass","Obs Mass","Error/ppm","Intensity")
        self.calibration_box.column("#0",width=0)
        self.calibration_box.column("Label",width=80)
        self.calibration_box.column("Ref Mass",width=80)
        self.calibration_box.column("Obs Mass",width=80)
        self.calibration_box.column("Error/ppm",width=80)
        self.calibration_box.column("Intensity",width=80)
        self.calibration_box.heading("#0",text="")
        self.calibration_box.heading("Label",text="Label")
        self.calibration_box.heading("Ref Mass",text="Ref Mass")
        self.calibration_box.heading("Obs Mass",text="Obs Mass")
        self.calibration_box.heading("Error/ppm",text="Error/ppm")
        self.calibration_box.heading("Intensity",text="Intensity")
        self.calibration_box.insert(parent="",index=0,iid=0,text="",values=("select","a","reference","mass","list"))
        self.calibration_box.grid(row=2,column=0,columnspan=4,rowspan=6)
        self.removeselected_button=tk.Button(self.callibration_tab,text="remove",command=self.remove_calibration_point).grid(row=2,column=4)
        self.calibration_type_label=tk.Label(self.callibration_tab,text="calibration type").grid(row=9,column=0)
        self.calibration_type=ttk.Combobox(self.callibration_tab,values=["Ledford","Quadratic"])
        self.calibration_type.grid(row=9,column=1)
        self.clibration_accept_button=tk.Button(self.callibration_tab,text="accept",command=self.accept_cal).grid(row=3,column=4)
        self.calibration_type.set("ledford")
        
        
    def accept_cal(self):
        for i in self.d1masslistbox.get_children():
            self.d1masslistbox.delete(i)
        for i in range(len(self.d1peaks)):
            self.d1peaks[i][0]=quadratic(self.d1peaks[i][0],*self.calparams)
            try:
                self.d1masslistbox.insert("","end",str(self.peaksd1[i][1]),values=[self.peaksd1[i][0],self.peaksd1[i][1]])
            except:
                o=2
        for i in range(len(self.fragmentaxis)):
            self.fragmentaxis[i]=quadratic(self.fragmentaxis[i],*self.calparams)
        
        self.callibration_tab.destroy()

    def remove_calibration_point(self):
        choices=[]
        focussed=self.calibration_box.focus()
        self.removed.append(float(self.calibration_box.item(focussed)["values"][1]))
        for i in range(len(self.choice_of_peaks)):
            if self.choice_of_peaks[i][1] not in self.removed:
                choices.append(self.choice_of_peaks[i])
        self.update_box(choices)
        
    def select_cal_list(self):
        self.reflistfile="reflists/"+self.callibration_combobox.get()
        self.refs=import_ref_file(self.reflistfile)
        cal=1
        if cal==0:
            print("something went wrong")
        else:
            self.choice_of_peaks=self.getcalchoices()            
            self.update_box(self.choice_of_peaks)
            
    def select_cal_list_with_variable(self,reflist):
        self.reflistfile="reflists/"+reflist
        self.refs=import_ref_file(self.reflistfile)
        if not self.spectrum:
            focussed=self.calibration_box.focus()
            self.calibration_box.item(0,values=("select","a","Mass","Spec","File"))
        else:
            self.choice_of_peaks=self.getcalchoices()            
            self.update_box(self.choice_of_peaks)    
            
    def getcalchoices(self):
        refmasses=np.transpose(self.refs)[1]
        #self.d1peaks=peak_picker_guassian1d([self.linex[10:-10],self.liney[10:-10]])
        calmasses=np.transpose(self.d1peaks)[0]
        labels=np.transpose(self.refs)[0]
        calchoices=[]
        for mass in  range(len(refmasses)):
            closesti=np.argmin(np.abs(np.subtract(calmasses,float(refmasses[mass]))))
            if abs(mass_error(float(mass), calmasses[closesti]))<100:
                calchoices.append([labels[mass],float(refmasses[mass]),calmasses[closesti],mass_error(float(refmasses[mass]), float(calmasses[closesti])),self.d1peaks[closesti][1]])
        return calchoices
        
    def update_box(self,entries):
        xdata=[]
        ydata=[]
        for i in self.calibration_box.get_children():
            self.calibration_box.delete(i)
        v=(self.calibration_type.get())
        if self.calibration_type.get()=="ledford":
            entriest=np.transpose(entries)
            xdata,ydata=entriest[1],entriest[0]
            params,params_covarience=optimize.curve_fit(quadratic,xdata,ydata)#,p0=[1,1,1])
            self.calparams=params
        elif self.calibration_type.get()=="Quadratic":
            for i in range(len(entries)):
                        xdata.append(entries[i][2])
                        ydata.append(entries[i][1])
            ML1,ML2,ML3,=float(self.params["ML1"]),float(self.params["ML2"]),float(self.params["ML3"])
            params,params_covarience=optimize.curve_fit(bruker_quadratic,xdata,ydata,p0=[ML1,ML2,ML3])
            self.calparams=params
        for i in range(len(entries)):
            a=i
            if v=="ledford":
                self.calibration_box.insert(parent="",index=a,iid=a,values=(entries[a][0],entries[a][1],quadratic(entries[a][1],*params),mass_error(entries[a][1],quadratic(entries[a][0],*params)),entries[a][3]))
            elif v=="Quadratic":
                self.calibration_box.insert(parent="",index=a,iid=a,values=(entries[a][0],entries[a][1],bruker_quadratic(entries[a][1],*params),mass_error(entries[a][1],bruker_quadratic(entries[a][0],*params)),entries[a][3]))

def gaussian(x, amp, cent,width):
    return amp *np.exp(-((x-abs(cent))**2)/(2*width**2))

def fitgaussian(x,y,po):
    params,covar=scipy.optimize.curve_fit(gaussian,x,y,p0=po)
    return params

def COM_calculator(spectrum,position):
    COM=((spectrum[1][position-3]*spectrum[0][position-3])+(spectrum[1][position-2]*spectrum[0][position-2])+(spectrum[1][position-1]*spectrum[0][position-1])+(spectrum[1][position+1]*spectrum[0][position+1])+(spectrum[1][position+2]*spectrum[0][position+2])+(spectrum[1][position+3]*spectrum[0][position+3])+(spectrum[1][position]*spectrum[0][position]))/(spectrum[1][position-1]+spectrum[1][position-2]+spectrum[1][position-3]+spectrum[1][position]+spectrum[1][position+1]+spectrum[1][position+2]+spectrum[1][position+3])
    return COM

def find_noise(y):
    for i in range(3):
        y=scipy.signal.savgol_filter(y, 3, 1)
    noise=np.std(y)*3*2.5
    return noise
def find_noise2(spectrum):
    #spectrum=spectrum.flatten("C")
    spectrum=np.sort(spectrum,axis=None)
    noise=np.mean(spectrum[:int((7*len(spectrum))/8)])*5
    return noise
def mass_error(theomass,obsmass):
    #print(obsmass)
    me=((float(obsmass)-float(theomass))/float(theomass))*1000000
    return me
    
def gaussian(x, amp, cent,width):
    return amp *np.exp(-((x-cent)**2)/(2*width**2))
def peak_picker_guassian(spectrum):
    peaks=[]
    windowedge=3
    noise=find_noise2(spectrum[2])*5   
    above_noise_positions=np.where(spectrum[2]>noise) # gives results as coords [y,x]
    for i in range(len(above_noise_positions[0])):
        if above_noise_positions[0][i]>windowedge+10 and above_noise_positions[1][i]>windowedge+10 and above_noise_positions[0][i]<len(spectrum[1])-windowedge-1 and above_noise_positions[1][i]<len(spectrum[0])-windowedge-1:
            if spectrum[1][above_noise_positions[0][i]]>200 and spectrum[1][above_noise_positions[0][i]]<3000 and spectrum[0][above_noise_positions[1][i]]>200 and spectrum[0][above_noise_positions[1][i]]<3000: 
                x,y=above_noise_positions[1][i],above_noise_positions[0][i]
                inten=spectrum[2][y,x]
                if inten==np.max(spectrum[2][y-10:y+10,x-10:x+10]):
                    try:
                        xpeak=list(fitgaussian(spectrum[0][x-10:x+10],spectrum[2][y,x-10:x+10],[spectrum[2][y,x],spectrum[0][x],0.2]))[1]
                        ypeak=list(fitgaussian(spectrum[1][y-10:y+10],spectrum[2][y-10:y+10,x],[spectrum[2][y,x],spectrum[1][y],0.2]))[1]
                        peaks.append([xpeak,ypeak,spectrum[2][above_noise_positions[0][i]][above_noise_positions[1][i]]])
                    except:
                        j=1
    #print(len(peaks))
    return peaks

def contourdataretrieval(spectrum):
    x=[]
    y=[]
    windowedge=3
    noise=find_noise(spectrum)*5
    print(noise)
    above_noise_positions=np.where(spectrum>noise)
    print(len(above_noise_positions[1]))
    for i in range(len(above_noise_positions[0])):
        if above_noise_positions[0][i]>windowedge+20 and above_noise_positions[1][i]>windowedge+20 and above_noise_positions[0][i]<len(spectrum[1])-windowedge-1 and above_noise_positions[1][i]<len(spectrum[0])-windowedge-1:
            if spectrum[above_noise_positions[0][i]][above_noise_positions[1][i]] == np.max(spectrum[above_noise_positions[0][i]-windowedge:above_noise_positions[0][i]+windowedge+1,above_noise_positions[1][i]-windowedge:above_noise_positions[1][i]+windowedge+1]):
               y.append(above_noise_positions[0][i])
               x.append(above_noise_positions[1][i])
    return x,y
                
    
def peak_picker_guassian1d(spectrum):
    peaks=[]
    windowedge=3
    noise=find_noise2(spectrum[1])*3
    above_noise_positions=np.where(spectrum[1]>noise) # gives results as coords [y,x]
    for i in range(len(above_noise_positions[0])):
        if above_noise_positions[0][i]>windowedge+10 and above_noise_positions[0][i]<len(spectrum[0])-windowedge-1:
            if spectrum[0][above_noise_positions[0][i]]<3000 and spectrum[0][above_noise_positions[0][i]]>200:
                x=above_noise_positions[0][i]
                #print(x)
                inten=spectrum[1][x]
                if inten==np.max(spectrum[1][x-10:x+10]):
                    if inten>((np.sum(spectrum[1][x-10:x+10])-inten)/19)+(np.std(spectrum[1][x-10:x-1])*3):
                        try:
                            xpeak=fitgaussian(spectrum[0][x-5:x+5], spectrum[1][x-5:x+5],po=[spectrum[1][x],spectrum[0][x],0.02])[1]
                            peaks.append([float(abs(xpeak)),spectrum[1][x]])
                        except:
                            print("couldnt work it out for:",spectrum[0][x])
    return peaks

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
def peak_picker_basic_3d_optimised(spectrum):
    peaks=[]
    windowedge=3
    noise=find_noise(spectrum[2])*5
    above_noise_positions=np.where(spectrum[2]>noise)
    print(len(above_noise_positions[1]))
    for i in range(len(above_noise_positions[0])):
        if above_noise_positions[0][i]>windowedge+20 and above_noise_positions[1][i]>windowedge+20 and above_noise_positions[0][i]<len(spectrum[1])-windowedge-1 and above_noise_positions[1][i]<len(spectrum[0])-windowedge-1:
            if spectrum[2][above_noise_positions[0][i]][above_noise_positions[1][i]] == np.max(spectrum[2][above_noise_positions[0][i]-windowedge:above_noise_positions[0][i]+windowedge+1,above_noise_positions[1][i]-windowedge:above_noise_positions[1][i]+windowedge+1]):
                comx=np.divide(np.sum(np.multiply(spectrum[2][above_noise_positions[0][i],above_noise_positions[1][i]-windowedge:above_noise_positions[1][i]+windowedge+1],spectrum[0][above_noise_positions[1][i]-windowedge:above_noise_positions[1][i]+windowedge+1])),np.sum(spectrum[2][above_noise_positions[0][i],above_noise_positions[1][i]-windowedge:above_noise_positions[1][i]+windowedge+1]))
                comy=np.divide(np.sum(np.multiply(spectrum[2][above_noise_positions[0][i]-windowedge:above_noise_positions[0][i]+windowedge+1,above_noise_positions[1][i]],spectrum[1][above_noise_positions[0][i]-windowedge:above_noise_positions[0][i]+windowedge+1])),np.sum(spectrum[2][above_noise_positions[0][i]-windowedge:above_noise_positions[0][i]+windowedge+1,above_noise_positions[1][i]]))
                peaks.append([comx,comy,spectrum[2][above_noise_positions[0][i]][above_noise_positions[1][i]]])
    print(len(peaks))
    return peaks

#calibration functions
def quadratic(x,A,B,C):
   return (B*(x**2))+(A*x)+C

def find_closest_indices(A, B):
    # Reshape arrays to make broadcasting work
    A=np.array(A)
    B=np.array(B)
    A_reshaped = A.reshape(-1, 1)
    B_reshaped = B.reshape(-1, 1)   
    differences = cdist(A_reshaped, B_reshaped) # Calculate absolute differences between A and B using cdist
    closest_indices = np.argmin(differences, axis=0)# Find indices with minimum differences for each point in B

    return closest_indices


def bruker_quadratic(x,ML1,ML2,ML3):
    return ((-ML1+np.sqrt((ML1**2-(4*ML3*(ML2-x)))))/(2*(ML2-x)))

def findkmdarray(refmass,miarray):
    kmds=[]
    mz=miarray[0]
    mkmd=math.floor(refmass)/refmass
    for i in mz:
        mkmr=i*mkmd
        kmds.append(math.ceil(mkmr)-mkmr)
    return kmds

def import_ref_file(ref_file):
    with open(ref_file,"r") as f:
        ref_table=[]
        ref_file_reader=csv.reader(f,delimiter="\t")
        for i in ref_file_reader:
            x=list(i)
            print(i)
            ref_table.append([x[0],float(x[1]),x[2]])
    return ref_table
if __name__=="__main__":
    test=showinggui()