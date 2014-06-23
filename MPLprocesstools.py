# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 10:54:41 2014

@author: dashamstyr
"""

def buffered_array(data,(x,y)):
    import numpy as np
    #create buffer around dataframe
    (rows,columns) = np.shape(data)
    
    b = int(np.ceil(y/2))
    t = y-b
    #create buffered array for calculating mean and std
    if x > 1:
        l = int(np.ceil(x/2))
        r = x-l
        newsize = (rows+x,columns+y)
        newarray = np.empty(newsize)
        (newrows,newcolums) = newarray.shape
        newarray[:l,b:-t] = data[0,:]
        newarray[l:-r,b:-t] = data
        newarray[-r:,b:-t] = data[-1,:]
    else:
        l=0
        r=0
        newsize = (rows,columns+y)
        newarray = np.empty(newsize)
        (newrows,newcolums) = newarray.shape
        newarray[:,b:-t] = data

    newarray[:,:b] = newarray[:,b:2*b]
    newarray[:,-t:] = newarray[:,-2*t:-t] 
    
    return newarray

def SNR_mask_depol(mplin,**kwargs):
    from copy import deepcopy
    
    SNRthreshold=kwargs.get('SNRthreshold',3)
    numprofs=kwargs.get('numprofs',1)
    bg_alt=kwargs.get('bg_alt',[])
    nopassval=kwargs.get('nopassval',float('nan'))
    inplace=kwargs.get('inplace',False)
    recalc=kwargs.get('recalc',False)
    datatype=kwargs.get('datatype','data')
    
    if recalc or not mplin.SNR:
        mplin = mplin.calculate_SNR(bg_alt,numprofs,datatype=['data'])
  
    #start by creating mask where areas that fall below SNRthreshold are zeroed out
    SNRmask = mplin.SNR[datatype][0]>=SNRthreshold
    
    if inplace:
        mplout=mplin
    else:
        mplout=deepcopy(mplin)
           
    if mplout.depolrat:
        mplout.depolrat[0]=mplin.depolrat[0]*SNRmask
        mplout.depolrat[0].replace(0,nopassval,inplace=True)
    return mplout 
    
def SNR_mask_all(mplin,**kwargs):
    from copy import deepcopy
    
    SNRthreshold=kwargs.get('SNRthreshold',3)
    numprofs=kwargs.get('numprofs',1)
    bg_alt=kwargs.get('bg_alt',[])
    nopassval=kwargs.get('nopassval',float('nan'))
    inplace=kwargs.get('inplace',False)
    recalc=kwargs.get('recalc',False)
    datatype=kwargs.get('datatype','data')
    if recalc or not mplin.SNR:
        mplin = mplin.calculate_SNR(bg_alt,numprofs,datatype=['data'])
  
    #start by creating mask where areas that fall below SNRthreshold are zeroed out
    SNRmask = mplin.SNR[datatype][0]>=SNRthreshold
    
    if inplace:
        mplout=mplin
    else:
        mplout=deepcopy(mplin)
    
    for n in range(mplout.header['numchans'][0]):
        mplout.data[n]=mplin.data[n]*SNRmask
        mplout.data[n].replace(0,nopassval,inplace=True)
        if mplout.rsq:
            mplout.rsq[n]=mplin.rsq[n]*SNRmask
            mplout.rsq[n].replace(0,nopassval,inplace=True)        
        if mplout.NRB:
            mplout.NRB[n]=mplin.NRB[n]*SNRmask
            mplout.NRB[n].replace(0,nopassval,inplace=True)        
    if mplout.depolrat:
        mplout.depolrat[0]=mplin.depolrat[0]*SNRmask
        mplout.depolrat[0].replace(0,nopassval,inplace=True)
    return mplout  
    
def NRB_mask_create(dfin,**kwargs):
    
    """
    generates threshold altitudes to avoids spurious results by removing all 
    data beyond strong signal spikes
    
    """
    import numpy as np
    import pandas as pan

    NRBthreshold=kwargs.get('NRBthreshold',3)
    NRBmin=kwargs.get('NRBmin',0.5)
    minalt=kwargs.get('minalt',150)
    numprofs=kwargs.get('numprofs',1)
    winsize=kwargs.get('winsize',5)
    
    #start by creating array of threshold altitudes and masking NRB copol
    #creates a new array with buffers to account for numprofs, winsize
    
    data = dfin.values
    altrange=dfin.columns.values
    (rows,columns) = data.shape
    minalt_index=np.where(altrange>=minalt)[0][0]
    newarray = buffered_array(data,(numprofs,winsize))
    (newrows,newcolums) = newarray.shape

    #set default values for cutoff to maximum altitude 
    threshalts=np.ones(len(dfin.index))*altrange[-1]
   
    for r in range(rows):
        tempprof=np.mean(newarray[r:r+numprofs],axis=0)
        for c in np.arange(minalt_index,columns):
            tempval = np.mean(tempprof[c:c+winsize])
            if tempval >= NRBthreshold:
                for c2 in np.arange(c,columns):
                    tempval = np.mean(tempprof[c2:c2+winsize])
                    if tempval <= NRBmin:
                        threshalts[r]=altrange[c2]
                        break 
    threshseries=pan.Series(data=threshalts,index=dfin.index)
    return threshseries

def NRB_mask_apply(dfin,threshseries,nopassval=float('nan'),inplace=True):
    
    if inplace:
        dfout=dfin
    else:
        from copy import deepcopy
        dfout=deepcopy(dfin)
    altvals = dfin.columns.values    
    for r in dfin.index:
        tempval=[x for x in altvals if x>=threshseries.ix[r]][0]
        dfout.ix[r,tempval:]=nopassval
    return dfout

def NRB_mask_all(MPLin,**kwargs):
    """
        uses a list of threshold altitudes, or generates one based on kwargs
        and applies it to all data sets within an MPL class object
    """
    import numpy as np
    
    threshseries=kwargs.get('threshseries',[])
    NRBthreshold=kwargs.get('NRBthreshold',3)
    NRBmin=kwargs.get('NRBmin',0.5)
    minalt=kwargs.get('minalt',150)
    numprofs=kwargs.get('numprofs',1)
    winsize=kwargs.get('winsize',5)
    nopassval=kwargs.get('nopassval',np.nan)
    inplace=kwargs.get('inplace',True)
    
    if inplace:
        MPLout=MPLin
    else:
        MPLout=MPLin.copy()
    
    if not threshseries:
        threshkwargs= {'NRBthreshold':NRBthreshold,'NRBmin':NRBmin,'minalt':minalt,
                       'numprofs':numprofs,'winsize':winsize,'nopassval':nopassval,
                       'inplace':inplace}
        try:
            threshseries=NRB_mask_create(MPLout.NRB[0],**threshkwargs)
        except IndexError:
            print "NRB Mask can only work for MPL class object where NRB has been calculated!"
            return MPLout
    
    for n in range(MPLout.header['numchans'][0]):
        MPLout.data[n]=NRB_mask_apply(MPLout.data[n],threshseries)
        
        if MPLout.rsq:
            MPLout.rsq[n]=NRB_mask_apply(MPLout.rsq[n],threshseries)
        if MPLout.NRB:
            MPLout.NRB[n]=NRB_mask_apply(MPLout.NRB[n],threshseries)
    
    if MPLout.depolrat:
        MPLout.depolrat[0]=NRB_mask_apply(MPLout.depolrat[0],threshseries)
    
    return MPLout
    
def slopecalc(prof, winsize = 5):
    import pandas as pan
    import numpy as np
    """
    Calculates slope of data for a single profile using a smoothing window of
    predetermined size
    
    inputs:
    prof:  a pandas series where index is altitude
    n:  number of consecutive values to use to calculate slope
    
    output:
    slope: output series,same size as input,with profile slopes
    """
    data = prof.values
    altrange = prof.index
    
    numvals = len(data)
    
    #calculate slopes of profile
      
    slopevals = np.empty_like(data)
    if winsize < 2:
        print "Sorry, need at least two values to calculate slope!"
        return
    elif winsize == 2:
        for i in range(numvals-1):
            slopevals[i] = (data[i+1]-data[i])/(altrange[i+1]-altrange[i])
        slopevals[-1] = slopevals[-2]
    else:
        #set up empty arrays for calculating slope with buffers on each end
        tempdat = np.empty(numvals+winsize)
        tempalt = np.empty_like(tempdat)
        #define buffer widths on left and right
        l = int(np.ceil(winsize/2))
        r = winsize-l
        #populate central segment with data and altrange
        tempdat[l:r] = data
        tempalt[l:r] = altrange
        #define window values for filling data and index buffers 
        lwinvals = data[:winsize]
        lwindex = altrange[:winsize]
        rwinvals = data[-winsize:]
        rwindex = altrange[-winsize:]
        #determine slope and intercept for left and right data windows      
        LA = np.array([lwindex,np.ones(winsize)])
        [lint,lslope]=np.linalg.lstsq(LA.T,lwinvals)
        RA = np.array([rwindex,np.ones(winsize)])
        [rint,rslope]=np.linalg.lstsq(RA.T,rwinvals)
        #calculate altitude index values for left and right buffers
        laltstep = altrange[1]-altrange[0]
        raltstep = altrange[-1]-altrange[-2]
        lold=altrange[0]
        rold=altrange[-1]
        for n in np.arange(l)[::-1]:
            tempalt[n]=lold-laltstep
            lold=tempalt[n]
        for m in np.arange(r):
            tempalt[m]=rold+raltstep
            rold=tempalt[m]
        #fill in data buffers
        tempdat[:l]=lint+lslope*tempalt[:l]
        tempdat[-r:]=rint+rslope*tempalt[-r:]
        
        #use linear regression to calculate slopes for sliding windows
        for v in range(numvals):
            windat = tempdat[v:v+winsize]
            winalt = tempalt[v:v+winsize]
            A = np.array([winalt,np.ones(winsize)])
            [inter,slopevals[v]]=np.linalg.lstsq(A.T,windat)
        
    slope = pan.Series(data=slopevals,index=altrange)
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
    import numpy as np
    
    BL_out = pan.Series(index = dfin.index)
    
    if algo=="slope":
        for i in dfin.index:
            tempslope = slopecalc(dfin.ix[i])
            tempslope = tempslope.ix[:maxalt]
            test = lambda s,k:s.ix[k] >= slope_thresh
            passfail = test(tempslope,tempslope.index)            
            slopepass = passfail.groupby(passfail).groups[True]
            tempindex = tempslope.ix[slopepass[0]:].index
            templist=[]
            for n in range(len(slopepass)):
                s=slopepass[n]
                if len(templist)==0:
                    indexcount=tempindex.get_loc(s)
                t=tempindex[indexcount]
                if s==t:
                    templist.append(s)
                    indexcount+=1
                    if len(templist)>=numvals:
                        BL_out.ix[i] = templist[0]
                        break
                else:
                    templist=[]
    
#    if algo=="value":
#        for i in dfin.index:
#            tempvals = dfin.ix[i].ix[:maxalt]
#            test = lambda s,k:s.ix[k] <= val_thresh
#            passfail = test(tempvals,tempvals.index)
#            valpass = passfail.groupby(passfail).groups[True]
#            tempindex = tempslope.ix[valpass[0]:].index
#            templist=[]
#            for v in range(len(valpass)):
#                v=valpass[n]
#                if len(templist)==0:
#                    indexcount=tempindex.get_loc(v)
#                t=tempindex[indexcount]
#                if v==t:
#                    templist.append(v)
#                    indexcount+=1
#                    if len(templist)>=numvals:
#                        BL_out.ix[i] = templist[0]
#                        break
#                else:
#                    templist=[]
    
    if algo=="combo":
        for i in dfin.index:
            tempslope = slopecalc(dfin.ix[i])
            tempslope = tempslope.ix[:maxalt]
            tempvals = dfin.ix[i].ix[:maxalt]
            slopetest = lambda s,k:s.ix[k] >= slope_thresh
            passfail_slope = slopetest(tempslope,tempslope.index)
            valtest = lambda s,k:s.ix[k] <= val_thresh
            passfail_vals = valtest(tempvals,tempvals.index)
            
            slopepass = passfail_slope.groupby(passfail_slope).groups[True]
            valpass = passfail_vals.groupby(passfail_vals).groups[True]            
            combopass= np.sort(list(set(slopepass+valpass)))
            
            tempindex = tempslope.ix[combopass[0]:].index
            templist=[]
            for c in range(len(combopass)):
                c=combopass[n]
                if len(templist)==0:
                    indexcount=tempindex.get_loc(c)
                t=tempindex[indexcount]
                if c==t:
                    templist.append(c)
                    indexcount+=1
                    if len(templist)>=numvals:
                        BL_out.ix[i] = templist[0]
                        break
                else:
                    templist=[]            
    
    BL_out.fillna(maxalt)
    
    return BL_out
    

#def SNR recursive

    # Step 1: Break df into chunks

    # Step 2: Generate SNR Mask
    # Step 3: Test chunks for SNR threshold
    # Step 4: increase chunk size for those that don't pass

#    
#    
#    
#
#def layer_detect(dfin, algo="slope",slope_thresh=[],val_thresh=[],numvals=1,altrange=[]):
#    """
#    Takes a pandas dataframe containing lidar profiles and uses a combination of
#    slopes and threshold levels to detect feature edges
#    
#    inputs:
#    dfin: a pandas dataframe wit datetime index and altitude columns
#    slope=[]: if slope is defined, this value is used as the threshold slope to demarcate layers
#    thresold=[] if 
#   """    


def find_layers(MPLin,**kwargs):
    """
    takes an MPL class object and process it, one profile at a time, to estimate
    bottom, peak,and top for each layer within the 2-D dataset
    
    inputs:
    MPLin = an MPL-class object to be proken inot layers
    
    kwargs:
    
    Outputs:
    layers = a three-column dataframe with date-time index, containing altitude values
    for bottom, peak,and top of layers
    """
    from scipy import signal
    import numpy as np
    import pandas as pan
    import process_tools as ptools
    import MPLtools as mtools
    
    #if MPLin does not have all necessary processed data,generate it
    wavelet=kwargs.get('wavelet',signal.ricker)
    widths=kwargs.get('widths',[2])
    layerwidth=kwargs.get('layerwidth',2)
    bg_alt=kwargs.get('bg_alt',[])
    noisethresh=kwargs.get('noisethresh',3)
    datatype=kwargs.get('datatype','data')

    MPLin=MPLin.calc_all()
    #use raw data from co-polarized channel (not r-squared corrected) to find layers
    rawdata=MPLin.data[0]
    SNR=MPLin.SNR[datatype][0]
    
    panelout=pan.Panel(major_axis=rawdata.index,minor_axis=['Base','Peak','Top'])
    
    for i in rawdata.index:
        tempprof=rawdata.ix[i]
        tempSNR=SNR.ix[i]
        #set baseline noise level based on 
        if bg_alt:
            tempsigma0=np.mean([s for s,ind in zip(tempSNR.values,tempSNR.index) if ind>=bg_alt])
        else:
            tempsigma0=np.mean(tempSNR.iloc[-100:])
        temp_cwt=signal.cwt(tempprof,wavelet,widths)
        tempmax=ptools.maxmin(temp_cwt,widths,np.greater)
        tempmin=ptools.maxmin(temp_cwt,widths,np.less)
        
        #use only width of 2 for now
        
        minloc=[minval[1] for minval in tempmin if minval[0]==layerwidth]
        maxloc=[maxval[1] for maxval in tempmax if maxval[0]==layerwidth]
        templayers=ptools.layer_filter(tempprof,maxloc,minloc,tempsigma0,noisethresh)
        for n in range(len(templayers)):
            panelname='Layer{0}'.format(n)
            panelout.loc[panelname,i]=templayers[n]
    
    return panelout
        
    
if __name__=='__main__':
    
    import MPLtools as mtools
    import MPLplot as mplot
    import numpy as np
    import pandas as pan
    
    os.chdir('C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Whistler-0330\Processed')
    
    mpltest = mtools.MPL()
    
    mpltest.fromHDF('201403300000-201403302300_proc.h5')
    
    
    p=find_layers(mpltest)   
    
    
    
#    copol_mask = SNR_mask_all(copol)
#    depol_mask=SNR_mask_all(depol)
#    
#    mpl.NRB[0]=copol_mask
#    mpl.depolrat[0]=depol_mask
#    kwargs = {'saveplot':False,'showplot':True,'verbose':True}
#    
#    mplot.doubleplot(mpl,**kwargs)
    
#    bl=boundary_layer_detect(copol,slope_thresh=-0.009,numvals=5)