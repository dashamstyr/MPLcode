# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 16:25:43 2014

@author: dashamstyr
"""
import os
import lidar_tools as ltools
import MPLtools as mtools
import MPLplot as mplot
import numpy as np
import pandas as pan
import lidar_tools as ltools
import MPLprocesstools as mproc

os.chdir('C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\Ucluelet Files\Processed')

MPLin = mtools.MPL()

MPLin.fromHDF('201405030000-201405030900_proc.h5')

wave=532.0
winsize=20 
varthresh=1

MPLin=MPLin.calc_all()    
NRBin=MPLin.NRB[0]
bg_alt=NRBin.columns.values[-100]
#Step 1: calculate molecular profile
z=NRBin.columns.values
Pmol=ltools.molprof(z,wave)
Pmol_cor=Pmol.ix['vals']
molmask=pan.DataFrame(index=NRBin.index,columns=NRBin.columns)


for proftime in NRBin.index:    
#Step 2:extract NRB profile
    tempprof=NRBin.ix[proftime]
    #Step 3: divide the two and average together until noise is tamped down
    temprat=Pmol_cor/tempprof
    bufferedprof=mproc.buffered_array(temprat.values,(1,winsize))
    coef=pan.Series(index=temprat.index)
    for n in range(len(temprat.values)):
        tempratvals=bufferedprof[n:n+winsize]
        coef.iloc[n]=np.mean(tempratvals)
    #Step 4: Result from 3 is profile of multiplying facto.Use this to calculate variance of residuals
#        rawvariance=(z**(-2.0)*(tempprof-(Pmol_cor/coef)))**2
    rawvariance=((1/(z/1000.0)**2.0)*(tempprof-(Pmol_cor/coef)))**2.0
    bufferedvariance=mproc.buffered_array(rawvariance.values,(1,winsize))
    variance=pan.Series(index=rawvariance.index)
    for n in range(len(rawvariance.values)):
        tempvarvals=bufferedvariance[n:n+winsize]
        variance.iloc[n]=np.mean(tempvarvals)
        
    
#        sigmacalcvals=[v for v,r in zip(tempdat.values,tempprof.index) if r>=bg_alt]
#        sigmacalcfilt=[x for x in sigmacalcvals if not np.isnan(x)]
#        sigmatemp=np.std(sigmacalcfilt)
    dataprof=MPLin.data[0].ix[proftime]
    SNR=mproc.calc_SNR(dataprof,bg_alt=bg_alt)
    sigmaprof=dataprof/SNR
    tempmask=[v<=varthresh*s**2 for v,s in zip(variance,sigmaprof)]
    molmask.ix[proftime]=tempmask



    
