# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 10:54:41 2014

@author: dashamstyr
"""

def calculate_SNR(dfin,boxsize = (10,10)):
    import pandas as pan
    import numpy as np
    """
    inputs:
    dfin = a pandas dataframe
    boxsize=(10,10)  a tuble defining the size of the box to calculate in (x,y)
    
    Takes in a pandas datafram and generates a dataframe of the same size containing
    2-D signal to noise ratio calculations.
    """
    
    #create buffer around dataframe
    data = dfin.values
    (rows,columns) = dfin.shape
    newsize = (rows+boxsize[0],columns+boxsize[1])
    newarray = np.empty(newsize)
    (newrows,newcolums) = newarray.shape
    
    l = int(np.ceil(boxsize[0]/2))
    r = boxsize[0]-leftbuffer
    t = int(np.ceil(boxsize[1]/2))
    b = boxsize[1]-topbuffer
    
    #create buffered array for calculating mean and std
    
    newarray[:l,t:-b] = data[:l,:]
    newarray[l:-r,t:-b] = data
    newarray[-r:,t:-b] = data[-r:,:]
    newarray[:,:t] = newarray[:,t:2*t]
    newarray[:,-b:] = newarray[:,-2*b:-b]
    
    #calculate SNR from mean and std
    SNRout = np.empty_like(data)
    for r in range(rows):
        for c in range(columns):
            window = newarray[r:r+boxsize[0],c:c+boxsize[1]]
            tempmean = np.mean(window)
            tempstd = np.std(window)
            SNRout[r,c] = tempmean/tempstd    
    
    dfout = pan.DataFrame(data = SNRout, index = dfin.index, columns = dfin.coulmns)
    
    return dfout

def calculate_slope(prof, n = 10):
    import pandas as pan
    import numpy as np
    """
    Calculates slope of data for a single profile using a smoothing window of
    predetermined size
    
    inputs:
    prof:  a pandas series where index is altitude
    n:  number of consecutive values to average
    
    output:
    slopeout: output series,same size as input,with profile slopes
    """
    data = prof.values
    altrange = prof.index
    
    numvals = len(data)
    
    #calculate slopes of profile
    rawslope = np.empty_like(data)
    for i in range(numvals-1):
        rawslope[i] = (data[i+1]-data[i])/(altrange[i+1]-altrange[i])
    rawslope[-1] = rawslope[-2]
    
    
    if n == 1:
        slope = pan.Series(data=rawslope,index=altrange)
    else:
        
        l = int(np.ceil(n/2))
        r = n-l
        tempslope = np.empty(numvals+n)
        tmepslope[:l] = rawslope[:l]
        tmepslope[l:-r] = rawslope
        tempslope[-r:] = rawslope[-r:]
        
        for i in range(numvals):
            smoothslope[i] = np.mean(tempslope[i:i+n])
        slope = pan.Series(data = smoothslope, index = altrange)
    
    return slope 

def boundary_layer_detect(dfin, algo="slope",slope_thresh=[],val_thresh=[],numvals=1,maxalt=2000):
    """
    Approximates the edge of the boundary layer using some combination
    of three algorithms:
    1) Negative slope threshold:  this defines the boundary layer as the lowest
    altitude the exceeds some negative slope threshold
    2) Value threshold:  this defines the top of the boundary layer as the lowest altitude
    for which the NRB dips below a value threshold
    3) Combination threshold:  Uses a combinaiton of the above two methods
    
    Inputs:
    dfin - A pandas dataframe containing a series of lidar profileswith altitude as the index
            and datetime as the column names
    algo - a string determining which algorithms to use.  Can be either:
            "slope" - purely slope threshold
            "value" - purely NRB threshold
            "combo" - combined algorithm
    slope_thresh - a floating point number defining the minimum slope to be used in the slope algorithm
    val_thresh - a floating point number defining the value threshold for the NRB value algorithm
    numvals - the number of consecutive values that must be below the slope or value threshold in order to be counted
    maxalt - an altitude in meters above which the algorithm is not longer applied
    
    Outputs:
    BL_out - a pandas series with datetime index and a boundary layer altitude in meters
    """
    import pandas as pan
    from itertools import groupby
    from operator import itemgetter
    
    
    BL_out = pan.Series(index = dfin.columns)
    
    if algo=="slope":
        for c in dfin.columns:
            tempslope = calculate_slope(dfin[c])
            tempslope = tempslope.ix[:maxalt]
            for k,g in groupby(enumerate(tempslope), lambda(i,s):s>=slope_thresh):
                temp = map(itemgetter(1),g)
                if len(temp) >= numvals:
                    BL_out[c] = temp[0]
                    break
    
    if algo=="value":
        for c in dfin.columns:
            tempvals = dfin[c].ix[:maxalt]
            for k,g in groupby(enumerate(tempvals), lambda(i,v):v<=val_thresh):
                temp = map(itemgetter(1),g)
                if len(temp) >= numvals:
                    BL_out[c] = temp[0]
                    break
    
    if algo=="combo":
        for c in dfin.columns:
            tempslope = calculate_slope(dfin[c])
            tempslope = tempslope.ix[:maxalt]
            tempvals = dfin[c].ix[:maxalt]
            tempdict = zip(tempslope,tempvals)
            for k,g in groupby(enumerate(tempslope), lambda(i,d):d[0]>=slope_thresh or d[1]<=val_thresh):
                temp = map(itemgetter(1),g)
                if len(temp) >= numvals:
                    BL_out[c] = temp[0]
                    break 
    
    BL_out.fillna(maxalt)
    
    return BL_out
    
def layer_detect(dfin, algo="slope",slope_thresh=[],val_thresh=[],numvals=1,altrange=[]):
    """
    Takes a pandas dataframe containing lidar profiles and uses a combination of
    slopes and threshold levels to detect feature edges
    
    inputs:
    dfin: a pandas dataframe wit datetime index and altitude columns
    slope=[]: if slope is defined, this value is used as the threshold slope to demarcate layers
    thresold=[] if 
    
    