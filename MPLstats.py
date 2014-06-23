# -*- coding: utf-8 -*-
"""
Created on Thu May 22 20:05:26 2014

@author: dashamstyr
"""
import os,sys,glob
import MPLtools as mtools
import MPLprocesstools as mproc
import MPLplot as mplot
import MPLbinit as binit
import MPL_hist2d as h2d
import time, datetime
import numpy as np

#find all .mpl or .h5 files between two dates

def filegetter(filedir=[],filetype='.mpl'):
    """
    takes a top-level directory and finds all relevant files 
    returns a list of filenames to iterate through
    
    inputs:
    filedir = top-level directory.  If emplty, interactive window is used
    filetype = file extension to look for (.mpl or .h5), defaults to .mpl
    """    
    
    if not filedir:
        filedir = mtools.set_dir('Sepect top-level folder to search')
    
    filelist=[]
    for root, _, filename in os.walk(filedir):
        if filename.endswith(filetype):
            filelist.append(os.path.join(root,filename))
    
    return filelist
        
def dataextractor(filename,**kwargs):
    """
    takes a file name, opens the file and extracts the MPL data into a dataframe
    and filters data based on kwargs returns an MPL class object
    
    inputs:
    filename = full-path from root to file to be extracted
    
    kwargs:
    altitudes = range of altitudes to sample
    starttime = datetime object showing earliest time to sample
    endtime = datetime object showing latest time to sample
    timestep = averaging time between profiles
    verbose = whether to print progress markers
    
    """
    
    altitudes = kwargs.get('altitudes',[])
    starttime = kwargs.get('starttime',[])
    endtime = kwargs.get('endtime',[])
    timestep = kwargs.get('starttime',[])
    verbose = kwargs.get('verbose',False)
    
    temp = mtools.MPL()
    if filename.endswith('.mpl'):        
        temp.fromMPL(filename)
    elif filename.endswith('.h5'):
        temp.fromHDF(filename)
    else:
        print "Sorry, cannot extract data from {0} files".format(filename.split('.')[-1])
    
    temptimes = temp.data[0].columns.values
    tempstart = temptimes[0]
    tempend = temptimes[-1]
    
    if starttime:
        if tempend < starttime:
            return
        else:
            newstart = starttime
    else:
        newstart = []
    
    if endtime:
        if tempstart > endtime:
            return
        else:
            newend = endtime
    else:
        newend = []
    
    mpl_out = temp.time_resample(starttime=newstart,endtime=newend,timestep=timestep,verbose=verbose)
    
    if altitudes:
        mpl_out = mpl_out.alt_resample(altitudes,verbose=verbose)
    
    return mpl_out        

def makebins(dfin,**kwargs):
    """
    takes a dataframe and generates a histogram of the data within based on kwargs
    returns a pandas series with index=bins and values=counts
    
    inputs
    ------------
    dfin = dataframe object to be histogrammed
    
    kwargs
    -----------
    minval = minimium value for histogram bins
    maxval = maximum value for histogram bins
    numbins = number of bins to use
    missinglowval = integer to assign to data below minval
    missinghighval = integer to assign to data above maxval
    
    """
    minval=kwargs.get('minval',dfin.min().min())  
    maxval=kwargs.get('maxval',dfin.max().max())
    numbins=kwargs.get('numbins',100)
    missinglowval=kwargs.get('missinglowval',-99999)
    missinghighval=kwargs.get('missinghighval',99999)
    
    histbins = binit.binit(minval,maxval,numbins,missinglowval,missinghighval)
    
    data_vec = np.ravel(dfin.values)
    histbins.do_bins(data_vec)
    
    return histbins

def multifile_hist(filedir=[],filetype='.mpl',histtype='copol',**kwargs):
    """
    opens all files within a directory, limited by kwargs, one by one, extracts
    data, and generates a histogram including data from all of them.  
    Values to be histogrammed and bin parameters are defined in kwargs
    
    inputs:
    
    filedir = top-level directory to search within (including sub-folders)
    filetype = file extension to search for ('.mpl' or '.h5')
    histtype = type of data to histogram (could be "copol","crosspol",
    "depolrat","SNR","layer_heights","layer_peaks")
    
    kwargs:
    
    altitudes = range of altitudes to scan
    starttime = earliest date/time to scan
    endtime = latest date/time to scan
    timestep = averaging time between profiles
    minval = minimum value for histograms (defaults to minimum of data set)
    maxval = maximum value for histograms (defaults to maximum of data set)
    numbins = number of bins to include (defaults to 100)
    missingLowVal = integer to assign to data below minval (defaults to -99999)
    missingHighVal = integer to assign to data above maxval (defaults to +99999)
    SNRmask = boolean to determine whether to apply SNR threshold mask
    SNRthresh = minimum SNR value to allow
    
    verbose = determines whether to print progress messages
    
    """
    altitudes=kwargs.get('altitudes',[])
    starttime=kwargs.get('starttime',[])
    endtime=kwargs.get('endtime',[])
    timestep=kwargs.get('timestep',[])
    minval=kwargs.get('minval',[])
    maxval=kwargs.get('maxval',[])
    numbins=kwargs.get('numbins',100)
    missingLowVal=kwargs.get('missingLowVal',-99999)
    missingHighVal=kwargs.get('missingHighVal',99999)
    SNRmask=kwargs.get('SNRmask',False)
    SNRthresh=kwargs.get('SNRthresh',1)
    SNRwin=kwargs.get('SNRwin',(5,5))
    verbose=kwargs.get('verbose',False)
    
    filelist=filegetter(filedir)
    binit_all=binit.binit(minval,maxval,numbins,missingHighVal,missingLowVal)   
    
    for f in filelist:
        mpl_temp = dataextractor(f,altitudes,starttime,endtime,timestep,verbose)        
        mpl_temp = mpl_temp.calc_all()        
        
        if histtype=="copol":
            df_temp=mpl_temp.NRB[0]
        elif histtype=="crosspol":
            df_temp=mpl_temp.NRB[1]
        elif histtype=="depolrat":
            df_temp=mpl_temp.depolrat[0]
        elif histtype=="SNR":
            df_temp=mpl_temp.SNR[0]
        
        if SNRmask:
            df_temp = mproc.SNR_mask_all(df_temp,SNRthresh,SNRwin)
        
        binit_temp = makebins(df_temp,minval,maxval,numbins,missingHighVal,missingLowVal)
        
        binit_all.append(binit_temp)
            
   return binit_all     
    
#filter by altitude, signal, depol, SNR
#choose stats type
#set up bins
#open each file and count in bins
#calculate stats mean, stdev, 
#combine with existing stats
#combine counts with existing counts
#

#