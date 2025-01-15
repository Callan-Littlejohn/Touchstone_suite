# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 10:28:43 2024

@author: FTICR_Kool_Kidz_PC
"""

import twoDMS_processing
import numpy

toprocess=["I:/15T/20240731_15T/2DMS_1M512K_nESI_Polyester_EID_0.275s_18.5bias_19.3lens_5scansDwell_000002_rfg3.d","I:/15T/20240731_15T/2DMS_1M513_nESI_Polyester_EID_0.275s_18.5bias_19.3lens_5scansDwell_000002.d"]
outdata="C:/Users/FTICR_Kool_Kidz_PC/Documents/testingtouch/batchstorage/"
outputpath=[]
for i in toprocess:
    name=i.split("/")[-1][:-2]
    outputpath.append(outdata+name)
for i in range(len(toprocess)):
    d2class=twoDMS_processing.d2spectrum(toprocess[i])
    d2class.process2d(urqrdrank=25)
    d2class.save2d(outputpath[i])