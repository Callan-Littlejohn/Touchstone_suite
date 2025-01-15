# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 13:08:02 2024

@author: FTICR_Kool_Kidz_PC
"""

import numpy as np
import matplotlib.pyplot as plot
from pyximport import install; install()
import array
import curqrd
sampling_rate=2**20
duration=0.1
no_zerofills=1

# file="C:/Users/FTICR_Kool_Kidz_PC/Documents/testingtouch/testingurqrd.metal/spectrum.npy"
# data=np.load(file)
# print(data)
def import_data(filename,n_scans,s_size): # taken directly from spike
    data=[]
    with open(filename,"rb") as f:
        for i in range(n_scans):
            indiv_scan=f.read(4*s_size)
            indiv_scan=array.array("l",indiv_scan)
            data.append(indiv_scan)
    return data
def gen_signal(freq,samplrate,duration, amplitude):
    x=np.linspace(0.0,duration,int(duration*samplrate),endpoint=False)
    freqs=x*freq
    y=amplitude*np.sin((2*np.pi)*freqs)
    return x,y

#data=import_data("H:/solarix 12T/20220107/20220107_HPmix_sod_000005.d/fid",1,2**23)[0]
#data=np.array(data)
a3,khz5=gen_signal(5000,sampling_rate,duration,1)
khz5p1=khz5+1
a2,khz4=gen_signal(4000,sampling_rate,duration,1)
a1,khz3=gen_signal(3000,sampling_rate,duration,1)
data=khz5+khz3+khz4

#plot.plot(data[0,1000: ])
#infile=