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

import numpy as np
import array
import os
from xml.dom import minidom
import time
import math
import matplotlib.pyplot as plot
import configparser
import matplotlib
import scipy.signal

#universal constants
da=1.66054e-27
echar=1.60217662e-19
emass=0.0005485833
font = {'family' : 'Helvetica',
        'weight' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)

class d2spectrum(): # this is the main class used to perform 2DMS processing and hold data
    def __init__(self,filename):
        self.filename=filename
        paramfile=[ name for name in os.listdir(filename) if os.path.isdir(os.path.join(filename, name)) ][0] # this kind of data should only have one sub folder so this is finding that subfolder which can be named many different things
        paramfile=filename+"\\"+paramfile+"\\apexAcquisition.method" #get method file
        self.params=read_param(paramfile) # get parameters
    
    def process2d(self,L_20=0,no_zerofillsx=1 ,no_zerofillsy=1,apodisation="Kaiser",beautify=False):
        #importing data
        beautification_array=[]
        t0=time.time() # 0 time
        if L_20==0: #default
            L_20=int(self.params["L_20"]) #kind of hacky way of allowing easier cutting
        serpath=os.path.join(self.filename,"ser") #the transients are located in the ser file
        data=np.asarray(import_data(serpath,L_20,int(self.params["TD"])),dtype=np.longdouble) # keeping as an np array makes life much easier (no need to transpose)
        t1=time.time() #import done time
        print("import done:",t1-t0)
        
        #apodisation generation
        if apodisation=="Kaiser": #this is what spike py uses
            afuncx=np.kaiser(int(self.params["TD"]),5)
            afuncy=np.kaiser(int(L_20),5)
        elif apodisation=="None": # no apodisation
            afuncx=np.ones(int(self.params["TD"]))
            afuncy=np.ones(int(L_20))
        else: #failsafe
            print("apodisation not recognised defaulting to none")
            afunc=1
        
        for i in range(len(data)): # process x dimension
            x=data[i]*afuncx # apodise .        operating on one slice at a time cuts down ram usage
            x=zerofilling(x,no_zerofillsx) # zerofill
            x=np.abs(np.fft.fft(x))#ft
            x=x[:len(data[i])]
            #x=np.sqrt(np.power(x.real,2)+np.power(x.imag,2))
            data[i]=np.abs(x) #feed back data to original space
        t2=time.time() #1st ft done time
        print("first ft complete:",t2-t0)
        for i in range(len(data[0])): #process y dimension
            x=data[:,i]*afuncy #apodise        operating on one slice at a time cuts down ram usage, operating on columns rather than transposing cuts down time and ram
            x=zerofilling(x,no_zerofillsy) # zerofill
            x=np.abs(np.fft.fft(x)) # ft
            x=x[:len(data[:,i])]
            if beautify==True:
                beautification_array.append(find_noise(x))
            #x=np.sqrt(np.power(x.real,2)+np.power(x.imag,2))
            data[:,i]=np.abs(x) #feed back data to original space
            if i%1E5==0:
                print(i) # really nice sanity check for operation
                #print(len(x[1:]))
        t3=time.time() # second ft done time
        print("second ft comlete:",t3-t0)
        if beautify==True:
            avenoise=np.mean(beautification_array)
            beautification_array=np.divide(beautification_array,avenoise)
            data=np.divide(np.power(np.divide(data,beautification_array),3),1E10)
            
            
        magnet=((float(self.params["MW_low"])*da)/echar)*2*np.pi*float(self.params["EXC_Freq_High"]) # determine magnetic field strength
        sampling_rate_x,transient_length_x,N_x=transient_params(float(self.params["MW_low"]),magnet,no_zerofillsx,int(self.params["TD"])) # get parameters for transient
        sampling_rate_y,transient_length_y,N_y=transient_params(ydimensionmasslow(float(self.params["IN_26"]),magnet) ,magnet,no_zerofillsy,L_20) # get parameters for transient
          # get x axis (m/z)
        x=np.fft.fftfreq(int(self.params["TD"])*(2**no_zerofillsx),1/sampling_rate_x)[:len(data[1])]
        y=np.fft.fftfreq(L_20*(2**no_zerofillsy),1/sampling_rate_y)[:len(data)]
        x=getmz(x[1:],magnet)
        y=np.add(y,float(self.params["EXC_Freq_Low"])) # frequency shift to account for bruker electronics
        y=getmz(y[1:],magnet) # get y axis (m/z)
        print(len(x),len(data[0]),len(y),len(data)) #sanity check
        self.spectrum=[x,y,data] #keep data in class for further processing
    
    def save2d(self,dirname,filetype="npy"):
        if filetype=="npy":
            self.save2dnp(dirname)
        elif filetype=="csv":
            self.save2dcsv(dirname)
    
    def save2dnp(self,dirname):
        dirname=dirname+".metal" #set filename with extension
        os.mkdir(dirname) #make the directory
        np.save(str(dirname+"/spectrum.npy"), self.spectrum[2])
        with open(str(dirname+"/axesinfo.ax"),"w") as f: #writing axes
            f.write("0,")
            f.write(",".join([str(self.spectrum[0][j]) for j in range(len(self.spectrum[0][1:]))]))
            f.write(",\n")
            f.write("0,")
            f.write(",".join([str(self.spectrum[1][j]) for j in range(len(self.spectrum[1]))]))
            f.close()
        configset=configparser.ConfigParser() 
        with open(str(dirname+"/parameters.p"),"w") as f: #setting up parameters
            configset["spectrum_information"]={"spectrum type":"2d","spectrum size": len(self.spectrum[2][0]),"ysize":len(self.spectrum[1])}
            configset["bruker_info"]={"TD":self.params["TD"],"ML1":self.params["ML1"],"ML2":self.params["ML2"],"ML3":self.params["ML3"],"SW_h":self.params["SW_h"]}
            configset.write(f) #write params ths is a good sanity check
        print("saved")
        
    def save2dcsv(self,dirname):
        dirname=dirname+".metal" #set filename with extension
        os.mkdir(dirname) #make the directory
        with open(str(dirname+"/spectrum.met"),"w") as f: #writing spectrum
            count=0
            for i in self.spectrum[2]:
                count=count+1
                x=",".join([str(j)for j in i]) #getting prepped for csv making
                f.write(x) #write line
                if count<len(self.spectrum[2]):            
                    f.write(",\n") # new line
            f.close() #close at end
        with open(str(dirname+"/axesinfo.ax"),"w") as f: #writing axes
            f.write("0,")
            f.write(",".join([str(self.spectrum[0][j]) for j in range(len(self.spectrum[0][1:]))]))
            f.write(",\n")
            f.write("0,")
            f.write(",".join([str(self.spectrum[1][j]) for j in range(len(self.spectrum[1]))]))
            f.close()
        configset=configparser.ConfigParser() 
        with open(str(dirname+"/parameters.p"),"w") as f: #setting up parameters
            configset["spectrum_information"]={"spectrum type":"2d","spectrum size": len(self.spectrum[2][0]),"ysize":len(self.spectrum[1])}
            configset["bruker_info"]={"TD":self.params["TD"],"ML1":self.params["ML1"],"ML2":self.params["ML2"],"ML3":self.params["ML3"],"SW_h":self.params["SW_h"]}
            configset.write(f) #write params ths is a good sanity check
        print("saved")
    
    def process_tims(self,start_v,end_v,steps,L_20=0,apodisation="Kaiser",no_zerofillsx=1):
        t0=time.time() # 0 time
        if L_20==0: #default
            L_20=int((abs(start_v)-abs(end_v))/steps)

        serpath=os.path.join(self.filename,"ser") #the transients are located in the ser file
        data=np.asarray(import_data(serpath,L_20,int(self.params["TD"]))) # keeping as an np array makes life much easier (no need to transpose)
        t1=time.time() #import done time
        print("import done:",t1-t0)
        
        #apodisation generation
        if apodisation=="Kaiser": #this is what spike py uses
            afuncx=np.kaiser(int(self.params["TD"]),5)
            afuncy=np.kaiser(int(L_20),5)
        elif apodisation=="None": # no apodisation
            afuncx=np.ones(int(self.params["TD"]))
            afuncy=np.ones(int(L_20))
        else: #failsafe
            print("apodisation not recognised defaulting to none")
            afunc=1
        for i in range(len(data)): # process x dimension
            x=data[i]*afuncx # apodise .        operating on one slice at a time cuts down ram usage
            x=zerofilling(x,no_zerofillsx) # zerofill
            x=np.abs(np.fft.fft(x))#ft
            x=x[:len(data[i])]
            data[i]=np.abs(x) #feed back data to original space
        t2=time.time() #1st ft done time
        print("first ft complete:",t2-t0)
        magnet=((float(self.params["MW_low"])*da)/echar)*2*np.pi*float(self.params["EXC_Freq_High"]) # determine magnetic field strength
        sampling_rate_x,transient_length_x,N_x=transient_params(float(self.params["MW_low"]),magnet,no_zerofillsx,int(self.params["TD"])) # get parameters for transient
        x=np.fft.rfftfreq(N_x,1/sampling_rate_x) # get x freq axis
        x=getmz(x[1:],magnet)
        y=np.linspace(start_v,end_v,L_20)
        self.spectrum=[x,y,data]
        

def import_data(filename,n_scans,s_size): # main function used to import data, adapted from spike-py
    data=[] #data array
    with open(filename,"rb") as f: # this data is stored as binary data
        for i in range(n_scans): #set number of scans
            indiv_scan=f.read(4*s_size) #get data
            indiv_scan=array.array("l",indiv_scan) #decode from little endians
            data.append(indiv_scan) #append to data array
    return data #send back data array

def read_param(method_name): # reading apexaqusiton xml
    xmldoc = minidom.parse(method_name)
    x = xmldoc.documentElement
    pp = {}
    children = x.childNodes
    for child in children:
        if (child.nodeName == 'paramlist'):
            params = child.childNodes
            for param in params:
                if (param.nodeName == 'param'):
                    k = str(param.getAttribute('name'))
                    for element in param.childNodes:
                       if element.nodeName == "value":
                           try:
                               v = str(element.firstChild.toxml())
                               #print v
                           except: 
                               v = "NC"
                    pp[k] = v
    return pp #returns dict of all parameters

def zerofilling(transient,no_zerofills): # zerofill
    for i in range(no_zerofills): #iterating allows for multiple zerofills
        transient=np.pad(transient,(0,len(transient)),"constant") # zerofill
    return transient #returns zerofilled transient

def find_noise(y):
    for i in range(3):
        y=scipy.signal.savgol_filter(y, 3, 1)
    noise=np.std(y)*3*2.5
    return noise

def transient_params(masslow,magnet,no_zerofills,s_size): 
    sampling_rate=2*getfreq(masslow,1,magnet) #work out nyqust rate and sampling rate
    #print(masslow,magnet)
    transient_length=s_size/sampling_rate # get transient length
    N=int(sampling_rate*(2**no_zerofills)*transient_length) #get number of data ponits
    return sampling_rate,transient_length,N #return all the params needed for transient processing

def getfreq(mass,charge,magnet): #get the frequency of a m/z
    mass=mass*da #work out charge in kg
    charge=charge*echar #work out charge in C
    freq=(magnet*charge)/(2*np.pi*mass) #determine frequency using base ICR equation . assumes no electric fields this is close then calibration eq get closer
    return freq # returns frequency

def getmz(freq,magnet): #get m/z of a freq
    mz=((magnet)/(freq*2*np.pi)) #get M/q in kg/C
    mz=mz*(echar/da) #convert to m/z
    return mz #return m/z

def ydimensionmasslow(in_26,magnet): #get the y dimension low from the in26 and magnet strength
    masslow= (magnet/(((1/in_26)/2)*2*math.pi))*(echar/da) #work out the mz
    return masslow #return y dimension mass low

def show_tims_heatmap(file):
    tim=d2spectrum(file)
    tim.process_tims(200, 50, 1)
    x=np.linspace(170,3000,3000-170+1)
    data=tim.spectrum
    data2=[]
    data2=(getbinnedspectrum(data[0], data[2], x))
    x=np.asarray(data[0])
    x=x.reshape(-1,2**9).mean(axis=1)
    l=[]
    for i in range(10):
        l.append((i+1)*1E6)
    for i in range(5):
        l.append((i+1)*2*1E7)
    
    plot.contourf(x,data[1],data2,levels=l[4:])
    plot.xlim(600,1100)
    plot.ylim(70,120)
    plot.xlabel("$\it{ m/z }$")
    plot.ylabel("Elution Voltage / V")
    plot.savefig("timscontourbsa.png",dpi=300)
    data2=[x,data[1],data2]
    return data2 

def getbinnedspectrum(mz,intensity,xaxis):
    m=2**9
    intensity=np.asarray(intensity)
    y=intensity.reshape(-1,m).sum(axis=1)
    y=y.reshape(len(intensity),int(len(intensity[0])/m))
    #mz=np.asarray(mz)
    #x=mz.reshape(-1,m).sum(axis=1)
    return y
