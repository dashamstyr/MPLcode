# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 10:54:41 2014

@author: dashamstyr
"""
import os,sys
import pandas as pan
import numpy as np
from copy import deepcopy
from itertools import groupby
from scipy import signal
import MPLtools as mtools
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import MPLplot as mplot
import operator
import inversiontools as itools
from scipy import optimize as opt
from copy import deepcopy

def calc_slope(prof, winsize = 10):
    
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
    altrange = np.asarray(prof.index.values,dtype='float')
    
    #Step 1: pad dataset to allow averaging
    
    leftpad = np.int(np.floor(winsize/2))
    rightpad = winsize-leftpad
      
    #Step 2: Calculate a linear fit to the data in the window
    
    slopes = np.empty(len(data)-winsize)
    for n in range(len(slopes)):       
        x = altrange[n:n+winsize]
        y = data[n:n+winsize]
        
        coeffs = np.polyfit(x,y,1,full=False)
        slopes[n] = coeffs[0]
        
    
    slopes = np.pad(slopes,(leftpad,rightpad),'edge')
    
    slope_out = pan.Series(slopes, index=altrange)
    
    
    return slope_out

def calc_SNR(prof,bg=[],bg_alt=[]):

    """
    inputs:
    prof = a pandas series
    bg = background signal level (stray light + dark current)
    bg_alt = altitude above which signal is assumed to be purely background
             if empty, topmost 100 data points are used1
    
    Calculates signal to noise ratios for mpl data
    """
        
    if not bg_alt:
        bg_alt=prof.index[-200]
    if not bg:
        bg = np.mean(prof.ix[bg_alt])
    
    SNRprof=pan.Series(np.empty_like(prof.values),index=prof.index)
    tempvals=[v for v,r in zip(prof.values,prof.index) if r>=bg_alt]
    tempfilt=[x for x in tempvals if not np.isnan(x)]
    sigmatemp=np.std(tempfilt)
    Ctemp=sigmatemp/np.mean(np.sqrt(np.abs(tempfilt)))
    SNR = lambda x: (x-bg)/(Ctemp*np.sqrt(np.abs(x)))
        
    SNRprof[:]=np.array([SNR(v) for v in prof.values]).clip(0)
        
    return SNRprof

def calc_sigma(prof,bg_alt=[],sigma0=[]):

    """
    inputs:
    prof = a pandas series
    bg = background signal level (stray light + dark current)
    bg_alt = altitude above which signal is assumed to be purely background
             if empty, topmost 100 data points are used1
    
    Calculates noise values for mpl data
    """
    
    if not sigma0:    
        if not bg_alt:
            bg_alt=prof.index[-200]
        
        tempvals=[v for v,r in zip(prof.values,prof.index) if r>=bg_alt]
        tempfilt=[x for x in tempvals if not np.isnan(x)]
        sigma0=np.std(tempfilt)
    Ctemp=sigma0/np.mean(np.sqrt(np.abs(tempfilt)))
    sigmacalc = lambda x: Ctemp*np.sqrt(np.abs(x))
        
    sigmaprof=np.array([sigmacalc(v) for v in prof.values]).clip(0)
        
    return sigmaprof

def calc_sigma_depolrat(copol,depol,bg_alt=[],sigma0=[]):
    #calculates std.dev of depol ratio using the concept of stdev of ratios of
    #dependent variables:  if R = A/B, (s(R)/R)^2 = (s(A)/A)^2+(s(B)/B)^2
    
    if not sigma0:
        if not bg_alt:
            bg_alt=copol.index[-200]
        
        copolsigma=calc_sigma(copol,bg_alt)
        depolsigma=calc_sigma(depol,bg_alt)
    else:
        copolsigma=calc_sigma(copol,sigma0)
        depolsigma=calc_sigma(depol,sigma0)
    
    depolrat=depol/copol
    depolratsigma=depolrat*np.sqrt((copolsigma/copol)**2+(depolsigma/depol)**2)
    
    return depolratsigma

def calc_SNR_depolrat(depolrat,**kwargs):
    
    copol=kwargs.get('copol',[])
    depol=kwargs.get('depol',[])
    depolratsigma=kwargs.get('depolratsigma',[])
    bg_alt=kwargs.get('bg_alt',[])
    sigma0=kwargs.get('sigma0',[])
    
    if not any(depolratsigma):
        if sigma0:
            depolratsigma=calc_sigma_depolrat(copol,depol,sigma0=sigma0)
        else:
            depolratsigma=calc_sigma_depolrat(copol,depol,bg_alt=bg_alt)
    
    SNR_out=depolrat/depolratsigma
    
    return SNR_out

def partdepolratcalc(depolin,beta_parallel,beta_mol,moldepolrat=0.0035):
    
    #default moldepolrat: narrow double filter allows only Cabannes line (see SPIE proc reference)
    A = moldepolrat/(1+moldepolrat)
    partdepolrat=(depolin*beta_parallel-A*beta_mol)/(beta_parallel-A*beta_mol)
    
    return partdepolrat

def buffered_array(data,(x,y)):

    #create buffer around dataframe
    datashape = np.shape(data)
    
    b = int(np.ceil(y/2))
    t = y-b
    
    if len(datashape)==1:
        rows=datashape[0]
        newsize=(rows+y)
        newarray=np.empty(newsize)
        (newrows)=newarray.shape
        newarray[b:-t]=data
        newarray[:b]=data[:b]
        newarray[-t:]=data[-t:]
        
    else:
        rows=datashape[0]
        columns=datashape[1]    
        #simply copy first and last values to fill in buffers
        if x > 1:
            l = int(np.ceil(x/2))
            r = x-l
            newsize=(rows+x,columns+y)
            newarray=np.empty(newsize)
            (newrows,newcolums)=newarray.shape
            newarray[:l,b:-t]=data[0,:]
            newarray[l:-r,b:-t]=data
            newarray[-r:,b:-t]=data[-1,:]
        else:
            l=0
            r=0
            newsize=(rows,columns+y)
            newarray=np.empty(newsize)
            (newrows,newcolums)=newarray.shape
            newarray[:,b:-t]=data
    
        newarray[:,:b]=newarray[:,b:2*b]
        newarray[:,-t:]=newarray[:,-2*t:-t] 
    
    return newarray

def SNR_mask_depol(mplin,**kwargs):
    
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
    
    if not any(threshseries):
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

def molecular_detect(MPLin,**kwargs):    
    
    wave=kwargs.get('wave',532.0)
    winsize=kwargs.get('winsize',30) 
    varthresh=kwargs.get('varthresh',1)
    savefile=kwargs.get('savefile',False)
    savefilename=kwargs.get('savefilename','testmolecular.h5')
    
    
    MPLin=MPLin.calc_all()    
    NRBin=MPLin.NRB[0]
    bg_alt=NRBin.columns.values[-200]
    #Step 1: calculate molecular profile
    z=NRBin.columns.values
    altstep=z[1]-z[0]  #assumes regular altitude steps throughout
    Pmol=itools.molprof(z,wave)
    Pmol_cor=Pmol['vals']
    panelout=pan.Panel(major_axis=NRBin.index,minor_axis=['Base','Top'])
    
    for proftime in NRBin.index:    
    #Step 2:extract NRB profile
        tempprof=NRBin.ix[proftime]
        #Step 3: divide the two and average together until noise is tamped down
        temprat=Pmol_cor/tempprof
        bufferedprof=buffered_array(temprat.values,(1,winsize))
        coef=pan.Series(index=temprat.index)
        for n in range(len(temprat.values)):
            tempratvals=bufferedprof[n:n+winsize]
            coef.iloc[n]=np.mean(tempratvals)
        #Step 4: Result from 3 is profile of multiplying facto.Use this to calculate variance of residuals
        rawvariance=((1/(z/1000.0)**2.0)*(tempprof-(Pmol_cor/coef)))**2.0
        bufferedvariance=buffered_array(rawvariance.values,(1,winsize))
        variance=pan.Series(index=rawvariance.index)
        for n in range(len(rawvariance.values)):
            tempvarvals=bufferedvariance[n:n+winsize]
            variance.iloc[n]=np.mean(tempvarvals)

        dataprof=MPLin.data[0].ix[proftime]
        sigmaprof=calc_sigma(dataprof,bg_alt=bg_alt)
        #Step 5: Regions where variance is below threshold multiple of noise varaince(sigma squared)
        #identified as molecular regions
        tempmask=pan.Series(z,index=[v<=varthresh*s**2 for v,s in zip(variance,sigmaprof)])
        tempgroups=tempmask.groupby(level=0)
        
        for g in tempgroups:
            if g[0]:
                tempalts=g[1]
                n=0
                for key,alt in groupby(enumerate(tempalts), lambda (i,x):i-(x-tempalts.iloc[0])/altstep):                    
                    layeralt=map(operator.itemgetter(1),alt)
                    panelname='Layer{0}'.format(n)
                    panelout.loc[panelname,proftime]=[layeralt[0],layeralt[-1]]
                    n+=1
    if savefile:
        store=pan.HDFStore(savefilename)
        store['molecular']=panelout
        
    return panelout
   
def molecular_detect_single(tempprof,dataprof,**kwargs):    
    
    wave=kwargs.get('wave',532.0)
    winsize=kwargs.get('winsize',30) 
    varthresh=kwargs.get('varthresh',1)
    savefile=kwargs.get('savefile',False)
    savefilename=kwargs.get('savefilename','testmolecular.h5')
    
    
    bg_alt=tempprof.index.values[-200]
    #Step 1: calculate molecular profile
    z=tempprof.index.values
    altstep=z[1]-z[0]  #assumes regular altitude steps throughout
    Pmol=itools.molprof(z,wave)
    Pmol_cor=Pmol['vals']
    
    #Step 3: divide the two and average together until noise is tamped down
    temprat=Pmol_cor/tempprof
    bufferedprof=buffered_array(temprat.values,(1,winsize))
    coef=pan.Series(index=temprat.index)
    for n in range(len(temprat.values)):
        tempratvals=bufferedprof[n:n+winsize]
        coef.iloc[n]=np.mean(tempratvals)
    #Step 4: Result from 3 is profile of multiplying facto.Use this to calculate variance of residuals
    rawvariance=((1/(z/1000.0)**2.0)*(tempprof-(Pmol_cor/coef)))**2.0
    bufferedvariance=buffered_array(rawvariance.values,(1,winsize))
    variance=pan.Series(index=rawvariance.index)
    for n in range(len(rawvariance.values)):
        tempvarvals=bufferedvariance[n:n+winsize]
        variance.iloc[n]=np.mean(tempvarvals)

    sigmaprof=calc_sigma(dataprof,bg_alt=bg_alt)
    #Step 5: Regions where variance is below threshold multiple of noise varaince(sigma squared)
    #identified as molecular regions
    tempmask=pan.Series(z,index=[v<=varthresh*s**2 for v,s in zip(variance,sigmaprof)])
    tempgroups=tempmask.groupby(level=0)
    dfout=pan.DataFrame(columns=['Base','Top'])
    basepoints=[]
    toppoints=[]
    for g in tempgroups:
        if g[0]:
            tempalts=g[1]
            n=0
            for key,alt in groupby(enumerate(tempalts), lambda (i,x):i-(x-tempalts.iloc[0])/altstep):                    
                layeralt=map(operator.itemgetter(1),alt)
                panelname='Layer{0}'.format(n)
                dfout.loc[panelname]=[layeralt[0],layeralt[-1]]
                basepoints.append((layeralt[0],tempprof.ix[layeralt[0]]))
                toppoints.append((layeralt[1],tempprof.ix[layeralt[1]]))
                n+=1
    
    
    fig=plt.figure()
    ax1=fig.add_subplot(2,1,1)
    ax1.plot(tempprof.index,tempprof.values)
    ax2=fig.add_subplot(2,1,2)
    ax2.plot(tempprof.index,varthresh*sigmaprof**2,tempprof.index,variance.values)
    ax2.set_ylim([0,max(varthresh*sigmaprof**2)*1.2])
    fig.canvas.draw()

def PBL_detect(MPLin,**kwargs):
    
    """
    """    
    wavelet=kwargs.get('wavelet',dog)
    widths=kwargs.get('widths',[4])
    layerwidth=kwargs.get('layerwidth',4)
    bg_alt=kwargs.get('bg_alt',[])
    noisethresh=kwargs.get('noisethresh',3)
    datatype=kwargs.get('datatype','NRB')
    layer_min=kwargs.get('layer_min',[])
    mol_min=kwargs.get('mol_min',[])

    MPLin=MPLin.calc_all()
    #typically use raw data from co-polarized channel (NRB) to find boundary layer
    if datatype=='data':
        rawdata=MPLin.data[0]
    elif datatype=='rsq':
        rawdata=MPLin.rsq[0]
    elif datatype=='NRB':
        rawdata=MPLin.NRB[0]
    elif datatype=='depolrat':
        rawdata=MPLin.depolrat[0]
    
    PBLout = pan.Series(index=rawdata.index)
    z=rawdata.columns.values
        
    for i in rawdata.index:
        tempprof=rawdata.ix[i]
        sigmaprof=calc_sigma(tempprof,bg_alt=bg_alt)
        if bg_alt:
            tempsigma0=np.mean(sigmaprof[bg_alt:])
        else:
            tempsigma0=np.mean(sigmaprof[-200:])
            
        tempcwt=signal.cwt(tempprof,wavelet,widths)
        tempmin=maxmin(tempcwt,widths,np.less)
        
        #use only width of layerwidth for now            
        minloc=[minval[1] for minval in tempmin if minval[0]==layerwidth]
        
        #PBL height can only be negative (min) values and edge effects are removed
        #by remving edges layerwidth in size
        
        layerwidthindex=np.where(widths==layerwidth)[0]
        CWTvals=tempcwt[layerwidthindex,:][0]
        CWTminvals=CWTvals[minloc]
        CWTminalts=z[minloc]
        
        tempmolht=mol_min.loc[i,'Base']
        
        try:
            templayerht=layer_min.loc[i,'Base']
        except AttributeError:
            templayerht=layer_min
            
        edgealt=z[layerwidth]
        
        if np.isnan(templayerht) or tempmolht<=templayerht:
            maxalt=tempmolht
        else:
            maxalt=templayerht
            
        try:
            PBLval=min([v[0] for v in zip(CWTminvals,CWTminalts) if edgealt<=v[1]<=maxalt])                
        except ValueError:
            PBLval=[]
            PBLout.ix[i]=maxalt
                
        if PBLval:
            PBLout.ix[i]=CWTminalts[np.where(CWTminvals==PBLval)]
            
    return PBLout

def maxmin(arrayin,widths,f):
    #find all local maxima or minina from each row of an array and return array
    #of index values for each max or min
    
    arrayout=[]
    for n in range(len(arrayin)):
        temp=signal.argrelextrema(arrayin[n],f)[0]        
        for t in temp:
            arrayout.append((widths[n],t))
        
    return arrayout

def dog(points,a):
    y_out=[]
    for x in range(points):
        x=x-points/2
        y=(np.exp(-x**2/(2*a**2))*-x)/(np.sqrt(2*np.pi)*a**3)
        y_out.append(y)
    return y_out    


def find_layers(MPLin,**kwargs):
    """
    takes an MPL class object and process it, one profile at a time, to estimate
    bottom, peak,and top for each layer within the 2-D dataset
    
    inputs:
    MPLin = an MPL-class object to be proken into layers
    
    kwargs:
    wavelet=type of wavelet to use to find layer edges.  default:signal.ricker
    widths=Range of wavelet widths to feed into CWT.  default:[2]
    layerwidth=wavelet width to use for final layer ID. default:2
    bg_alt=altitude to mark as particulate-free for background calc. defualt:[]
    noisethresh=threshold level to mark a layer above background noise. default:3
    cloudthresh=threshold signal level to mark a layer as cloud for (water,ice). default: (1.0,0.4)
    datatype=type of profile in MPL object to use for layers.  default:'data'
    savefile=boolean for whetehr to save results. default:False
    savefilename=name to save file under.  default:'testlayers.h5'
    
    
    Outputs:
    panelout = a pandas panel object with three axes:
        major-axis: datetime of individual profiles
        minor-axis: layer info ['Base','Peak','Top','Delta','Depol','Type',
                                'Sub-Type','Lidar_Ratio']
        columns: Layer number (e.g. 'Layer1')
    """
    
    #if MPLin does not have all necessary processed data,generate it
    wavelet=kwargs.get('wavelet',signal.ricker)
    widths=kwargs.get('widths',[2])
    CWTwidth=kwargs.get('CWTwidth',2)
    minwidth=kwargs.get('minwidth',4)
    bg_alt=kwargs.get('bg_alt',[])
    noisethresh=kwargs.get('noisethresh',3)
    cloudthresh=kwargs.get('cloudthresh',(1,0.4))
    datatype=kwargs.get('datatype','data')
    savefile=kwargs.get('savefile',False)
    savefilename=kwargs.get('savefilename','testlayers.h5')
    
    MPLin=MPLin.calc_all()
    #use raw data from co-polarized channel (not r-squared corrected) to find layers
    if datatype=='data':
        rawdata=MPLin.data[0]
        rawdepol=MPLin.data[1]
    elif datatype=='rsq':
        rawdata=MPLin.rsq[0]
        rawdepol=MPLin.rsq[1]
    elif datatype=='NRB':
        rawdata=MPLin.NRB[0]
        rawdepol=MPLin.NRB[1]
    
    rawdepolrat=MPLin.depolrat[0]
    
    panelout=pan.Panel(major_axis=rawdata.index,minor_axis=['Base','Peak','Top',
    'Delta','Depol','Type','Sub-Type','Lidar_Ratio'])
    
    for i in rawdata.index:
        tempprof=rawdata.ix[i]
        tempdepolprof=rawdepol.ix[i]
        tempdepolratprof=rawdepolrat.ix[i]
        z=tempprof.index
        #set baseline noise level based on 
        tempsigma0=np.mean(calc_sigma(tempprof,bg_alt))
        
        temp_cwt=signal.cwt(tempprof,wavelet,widths)
        tempmax=maxmin(temp_cwt,widths,np.greater)
        tempmin=maxmin(temp_cwt,widths,np.less)
        
        #use only width of 2 for now
        
        minloc=[minval[1] for minval in tempmin if minval[0]==CWTwidth]
        maxloc=[maxval[1] for maxval in tempmax if maxval[0]==CWTwidth]
        templayers=layer_filter(tempprof,tempdepolprof,maxloc,minloc,noisethresh,minwidth,sigma0=tempsigma0)
        
        for n in range(len(templayers)):
            indices=templayers[n]
            
            minalt=indices[0]
            peakalt=indices[1]
            maxalt=indices[2]
            delta=indices[3]
            meandepolrat=indices[4]
            panelname='Layer{0}'.format(n)
            layerdepolratprof=tempdepolratprof.ix[minalt:maxalt]
#            meandepolrat=np.mean(layerdepolratprof)
            peakval=tempprof.ix[peakalt]
            if peakval >= cloudthresh[0]:
                layertype='cloud'  
                layersubtype,layerratio=icewaterfilter(meandepolrat)
            elif peakval >= cloudthresh[1] and meandepolrat>0.25:
                layertype='cloud'  
                layersubtype,layerratio=icewaterfilter(meandepolrat)
            else:
                layertype='aerosol'
                layersubtype,layerratio=aerosoltypefilter(meandepolrat)
            panelout.loc[panelname,i,'Base']=minalt
            panelout.loc[panelname,i,'Peak']=peakalt
            panelout.loc[panelname,i,'Top']=maxalt
            panelout.loc[panelname,i,'Delta']=delta
            panelout.loc[panelname,i,'Depol']=meandepolrat
            panelout.loc[panelname,i,'Type']=layertype
            panelout.loc[panelname,i,'Sub-Type']=layersubtype
            panelout.loc[panelname,i,'Lidar_Ratio']=layerratio
    
    if savefile:
        store=pan.HDFStore(savefilename)
        store['layers']=panelout
        
    return panelout

def layer_filter(prof,depolprof,maxiloc,miniloc,thresh=3,minwidth=4,sigma0=[]):
    """
    takes a profile and a list of local maxima and minima from CWT analysis and calculates
    layer edges and peaks while filtering out peaks for which the delta from 
    edge to peak is less than some multiple of the shot noise from background and
    dark current
    
    once a layer is defined, it is then investigated for variations in depol ratio
    If significant variations exist, the layer is firther divided into sub-layers 
    based on these results
    
    inputs:
    prof - a pandas series represeting a single profile of lidar returns with altitude
    depolprof - a pandas series representing a profile of depol ratios with altitude
    maxiloc - a list of maximum index values from the CWT results at a given wavelet width
            represent the peaks of a given layer
    miniloc - a list of minimum index values from the CWT results at a given wavelet width
            represent the edges of a given layer
    sigma0 - baseline noise level for the profile, if empty it is calculated
    thresh - difference between peak and edge of a layer must exceeed this 
             multiple of sigma0 to be counted.  default: 3
    
    """
    #step 1:calculate noise floor, if not defined

    if not sigma0:
        sigma0=np.mean(calc_sigma(prof)[100:])
    #Step 2: Calculate profile values at each edge and peak
    layers=[]
    n=0
    nextminloc=0
    while n < len(maxiloc)-1:
        n+=1
        #note: due to buffering to obtain averages, first and last peaks are discarded
        peakloc=maxiloc[n]
        edge_below_list=[v for v in miniloc[nextminloc:] if v<peakloc]
        edge_above_list=[v for v in miniloc[nextminloc:] if v>peakloc]
        
        if not edge_above_list or not edge_below_list:
            continue
        #Step 3: Calculate delta signal between peak and lower edge (directly before)
        for edge_below in edge_below_list[::-1]:            
            delta_lower=prof.iloc[peakloc]-prof.iloc[edge_below]        
            #Step 4: Filter out false layers for which delta < thresh*signam0
            if delta_lower>thresh*sigma0:
                templowedge=edge_below
                break
            else:
                templowedge=[]
                #try to find upper edge where delta_upper exceeds threshold
        if templowedge:
            for edge_above in edge_above_list:
                delta_upper=prof.iloc[peakloc]-prof.iloc[edge_above]
                if delta_upper>thresh*sigma0:
                    #if upper edge is found, add indices of (lower,center,upper, maximum delta) to layers
#                    temppeakval=np.max(prof.iloc[templowedge:edge_above])
#                    temppeakloc=np.where(prof.values==temppeakval)
                    tempprof=prof.iloc[templowedge:edge_above]
                    if len(tempprof>=minwidth):
                        tempdepolprof=depolprof.iloc[templowedge:edge_above]
    #                    delta=max(delta_lower,delta_upper)
                        depol_layers=find_depollayers(tempprof,tempdepolprof,signalsigma0=sigma0)
                        
                        layers+=depol_layers
#                    layers.append((templowedge,temppeakloc,edge_above,max(delta_lower,delta_upper)))
                    try:
                        nextpeak=[p for p in maxiloc if p >edge_above][0]
                    except IndexError:
                        break
                    nextminloc=miniloc.index(edge_above)
                    n=maxiloc.index(nextpeak)-1
                    break
    return layers

def find_depollayers(copolprof,depolprof,**kwargs):
    widths=kwargs.get('widths',np.arange(2,5))
    layerwidth=kwargs.get('layerwidth',2)
    wavelet=kwargs.get('wavelet',signal.ricker)
    depolratsigma0=kwargs.get('depolsigma0',0.05)
    noisethresh=kwargs.get('noisethresh',1)
    signalsigma0=kwargs.get('signalsigma0',[])
#    prevlayernum=kwargs.get('prevlayernum',0)
    
    depolratMPL=depolprof/copolprof
    depolratprof=depolratMPL/(depolratMPL+1)
    depolratsigma=calc_sigma_depolrat(copolprof,depolprof,sigma0=signalsigma0)
    depolratSNR=calc_SNR_depolrat(depolratprof,depolratsigma=depolratsigma,signal0=signalsigma0)
    temp_cwt=signal.cwt(depolratSNR,wavelet,widths)
    tempmax=maxmin(temp_cwt,widths,np.greater)
    tempmin=maxmin(temp_cwt,widths,np.less)
    z=depolratprof.index
    
    minloc=[minval[1] for minval in tempmin if minval[0]==layerwidth]
    maxloc=[maxval[1] for maxval in tempmax if maxval[0]==layerwidth]
    
    edgeloc=np.sort(minloc+maxloc)
    edgealt=[z[v] for v in edgeloc]
    templayers=depol_filter(depolratprof,copolprof,edgealt,depolratsigma0,noisethresh)
    
    return templayers
    
#    panelout=pan.DataFrame(index=range(len(templayers)),columns=['Base','Peak','Top','Delta'])
#
#    for n in range(len(templayers)):
#        indices=templayers[n]
#        minalt=z[indices[0]]
#        maxalt=z[indices[1]]
#        delta=indices[2]
#        panelname='Layer{0}'.format(n+prevlayernum)
#        panelout.iloc[n,0]=minalt
#        panelout.iloc[n,1]=peakalt
#        panelout.iloc[n,2]=maxalt
#        panelout.iloc[n,3]=delta

def depol_filter(depolprof,signalprof,edgealts,sigma0,thresh=1):
    """
    
    """
    #Step 2: Calculate profile values at each edge and peak
    
    base=depolprof.index[0]
    n=0
    layers=[]
    while n<len(edgealts):
        edge=edgealts[n]
        try:
            nextedge=edgealts[n+1]
        except IndexError:
            nextedge=depolprof.index[-1]
           
        meanbelow=np.mean(depolprof.ix[base:edge])
        meanabove=np.mean(depolprof.ix[edge:nextedge])
        depoldelta=abs(meanabove-meanbelow)
        if depoldelta>=thresh*sigma0:
            if n==len(edgealts):
                bottom=edge
                top=nextedge
                meandepol=meanabove
            else:
                bottom=base
                top=edge
                meandepol=meanbelow
            temppeakval=np.max(signalprof.ix[bottom:top])
            temppeakloc=np.where(signalprof.values==temppeakval)
            temppeakalt=signalprof.index[temppeakloc].values[0]
            delta_below=signalprof.ix[temppeakalt]-signalprof.ix[bottom]
            delta_above=signalprof.ix[temppeakalt]-signalprof.ix[top]
            signaldelta=max(delta_below,delta_above)
            layers.append([bottom,temppeakalt,top,signaldelta,meandepol])
            base=edge
            n+=1
        else:
            n+=1
    if not layers:
        bottom=depolprof.index[0]
        top=depolprof.index[-1]
        meandepol=np.mean(depolprof)
        temppeakval=np.max(signalprof)
        temppeakloc=np.where(signalprof.values==temppeakval)
        temppeakalt=signalprof.index[temppeakloc].values[0]
        delta_below=signalprof.ix[temppeakalt]-signalprof.ix[bottom]
        delta_above=signalprof.ix[temppeakalt]-signalprof.ix[top]
        signaldelta=max(delta_below,delta_above)
        layers.append([bottom,temppeakalt,top,signaldelta,meandepol])
           
    return layers

def layerprofplot(profin,layersin,numlayer=30):
    
    z=profin.index
    vals=profin.values
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(vals,z)
    
    mcolors=['blue','red','green','yellow','orange','purple']
    if numlayer>len(layersin.ix['Base']):
        numlayer=len(layersin.ix['Base'])
        
    for n in range(numlayer):
        if n>(len(mcolors)-1):
            color=mcolors[int(np.floor(n/len(mcolors)))-1]
        else:
            color=mcolors[n]   
        templayer=layersin.iloc[:,n]
        if not np.isnan(templayer['Base']):
            ax.scatter(profin.ix[templayer.iloc[0]],templayer.iloc[0],c=color,marker='o')
            ax.scatter(profin.ix[templayer.iloc[1]],templayer.iloc[1],c=color,marker='x')
            ax.scatter(profin.ix[templayer.iloc[2]],templayer.iloc[2],c=color,marker='v')
            
    fig.canvas.draw()


def icewaterfilter(depolrat,**kwargs):
    waterthresh=kwargs.get('waterthresh',0.10)
    icethresh=kwargs.get('icethresh',0.25)
    
    if depolrat <= waterthresh:
        typeout='water'
        ratout=15.3
    elif depolrat <= icethresh:
        typeout='mixed'
        ratout= 25.0
    else:
        typeout='ice'
        ratout=50.0
    
    return typeout,ratout


def aerosoltypefilter(depolrat,**kwargs):
    smokethresh=kwargs.get('smokethresh',0.05)
    dustthresh=kwargs.get('dustthresh',0.15)
    
    if depolrat <= smokethresh:
        typeout='water_soluble'
        ratout=30.0
    elif depolrat <= dustthresh:
        typeout='smoke'
        ratout=70.0
    else:
        typeout='dust'
        ratout=40.0
    
    return typeout,ratout

    
def colormask_fromdict(mplin,pblin,molin,layersin):
    
    alts=mplin.NRB[0].columns
    times=mplin.NRB[0].index
    mask=pan.DataFrame(index=times,columns=alts)
    
    colordict={'molecular':0,
               'PBL':1,
               'ice':2,
               'water':3,
               'mixed':4,
               'dust':5,
               'smoke':6,
               'water_soluble':7,
               'unidentified':8}
    
    
    for t in times:
        tempprof=pan.Series(index=alts)
        
        for m in molin.items:
            tempmol=molin.ix[m,t]
            tempminalt=tempmol.ix['Base']
            tempmaxalt=tempmol.ix['Top']
            tempprof[(tempprof.index>=tempminalt) & (tempprof.index<=tempmaxalt)]=colordict['molecular']

        
        for l in layersin.items:
            templayer=layersin.ix[l,t]
            tempminalt=templayer.ix['Base']
            tempmaxalt=templayer.ix['Top']
            temptype=templayer.ix['Type']
            tempsubtype=templayer.ix['Sub-Type']
            try:
                tempprof[(tempprof.index>=tempminalt) & (tempprof.index<=tempmaxalt)]=colordict[tempsubtype]
            except KeyError:
                tempprof[(tempprof.index>=tempminalt) & (tempprof.index<=tempmaxalt)]=8
                
        PBLht=pblin.ix[t]
        tempprof[tempprof.index<=PBLht]=1

        tempprof.fillna(value=8,inplace=True)
        mask.ix[t]=tempprof
    
    return mask,colordict

def colormask_plot(maskin,colordict,**kwargs):
    #set color codes for different layers
    hours=kwargs.get('hours',['00','06','12','18'])
    fontsize=kwargs.get('fontsize',24)
    cbar_ticklocs=kwargs.get('cbarticklocs',np.arange(0,9)+0.5)
    altrange=kwargs.get('altrange',[])
    datetimerange=kwargs.get('datetimerange',[])
    SNRmask=kwargs.get('SNRmask',[])
    saveplot=kwargs.get('saveplot',True)
    plotfilepath=kwargs.get('plotfilepath',[])
    plotfilename=kwargs.get('plotfilename','testmaskfig.png')
    dpi = kwargs.get('dpi',100)

    cmapdict =  {'red':    ((0.0, 176.0/255.0, 176.0/255.0),
                            (0.105, 176.0/255.0, 255.0/255.0),
                            (0.21, 255.0/255.0, 255.0/255.0),
                            (0.33, 255.0/255.0, 0.0/255.0),
                            (0.44, 0.0/255.0, 186.0/255.0),
                            (0.55555, 186.0/255.0, 184.0/255.0),
                            (0.66666, 184.0/255.0, 0.0/255.0),
                            (0.77777, 0.0/255.0, 220.0/255.0),
                            (0.88888, 220.0/255.0, 192.0/255.0),
                            (1.0, 192.0/255.0, 192.0/255.0)),
        
                 'green':  ((0.0, 244.0/255.0, 244.0/255.0),
                            (0.105, 244.0/255.0, 69.0/255.0),
                            (0.21, 69.0/255.0, 255.0/255.0),
                            (0.33, 255.0/255.0, 0.0/255.0),
                            (0.44, 0.0/255.0, 85.0/255.0),
                            (0.55555, 85.0/255.0, 134.0/255.0),
                            (0.66666, 134.0/255.0, 0.0/255.0),
                            (0.77777, 0.0/255.0, 20.0/255.0),
                            (0.88888, 20.0/255.0, 192.0/255.0),
                            (1.0, 192.0/255.0, 192.0/255.0)),
        
                 'blue':   ((0.0, 230.0/255.0, 230.0/255.0),
                            (0.105, 230.0/255.0, 0.0/255.0),
                            (0.21, 0.0/255.0, 255.0/255.0),
                            (0.33, 255.0/255.0, 205.0/255.0),
                            (0.44, 205.0/255.0, 211.0/255.0),
                            (0.55555, 211.0/255.0, 11.0/255.0),
                            (0.66666, 11.0/255.0, 0.0/255.0),
                            (0.77777, 0.0/255.0, 60.0/255.0),
                            (0.88888, 60.0/255.0, 192.0/255.0),
                            (1.0, 192.0/255.0, 192.0/255.0))}    
                   
    maskmap=LinearSegmentedColormap('MaskMap',cmapdict)
    
    if altrange:
        maskin=maskin.loc[:,(maskin.columns>altrange[0]) & (maskin.columns<altrange[1])]
    
    if datetimerange:
        maskin=maskin[(maskin.index>datetimerange[0]) & (maskin.index<datetimerange[1])]
    
    datetime=maskin.index
    alts=maskin.columns
    
    if isinstance(SNRmask,pan.DataFrame):
        if altrange:
            SNRmask=SNRmask.loc[:,(SNRmask.columns>altrange[0]) & (SNRmask.columns<altrange[1])]
    
        if datetimerange:
            SNRmask=SNRmask[(SNRmask.index>datetimerange[0]) & (SNRmask.index<datetimerange[1])]
        
        maskin=maskin*SNRmask
        maskin.fillna(8,inplace=True)
    
    fig=plt.figure()
    ax1=plt.subplot2grid((37,60),(0,0),rowspan=30,colspan=60)
    cax=plt.subplot2grid((37,60),(30,0),rowspan=7,colspan=60)    
    image=ax1.imshow(maskin.T[::-1],cmap=maskmap,interpolation='none',vmin=0,vmax=9,aspect='auto')
    plt.tight_layout()
#    mplot.forceAspect(ax1,aspect=ar)
    mplot.dateticks(ax1, datetime, hours = hours,fsize=fontsize)
    ax1.set_xlabel('Hours [Local]',fontsize=fontsize+4)
    ax1.set_ylabel('Altitude [m]', fontsize=fontsize+4)
    mplot.altticks(ax1, alts[::-1], fsize = fontsize, tcolor = 'k')
#    divider = make_axes_locatable(ax)
#    cax = divider.append_axes("bottom", size="10%", pad=0.15)
    cbar1=fig.colorbar(image,cax=cax,orientation='horizontal')
#    cbar1.ax.set_autoscalex_on(False)
    cbar1.set_ticks(cbar_ticklocs)
    cbar1.ax.tick_params(bottom='off',top='off',labelsize=fontsize-4)
    
    sortedlabels=[s[0] for s in sorted(colordict.iteritems(), key=operator.itemgetter(1))]
    cbarlabels=cbar1.set_ticklabels(sortedlabels)
    
    if saveplot:
        if plotfilepath:
            if os.path.isdir(plotfilepath):
                savename=os.path.join(plotfilepath,plotfilename)
            else:
                os.mkdir(plotfilepath)
                savename=os.path.join(plotfilepath,plotfilename)
        fig.canvas.print_figure(savename,dpi = dpi, edgecolor = 'b', bbox_inches = 'tight') 
    fig.canvas.draw()
    
def findalllayers(filename,**kwargs):    
    timestep=kwargs.get('timestep','240S')
    maxalt=kwargs.get('maxalt',[])
    NRBmask=kwargs.get('NRBmask',True)
    NRBthresh=kwargs.get('NRBthresh',3)
    molthresh=kwargs.get('molthresh',1)
    molwinsize=kwargs.get('molwinsize',5)
    layernoisethresh=kwargs.get('layernoisethresh',0.4)
    cloudthresh=kwargs.get('cloudthresh',(1.0,0.20))
    layerdtype=kwargs.get('layerdtype','NRB')
    CWTwidth=kwargs.get('CWTwidth',2)
    layerminwidth=kwargs.get('minwidth',4)
    layerCWTrange=kwargs.get('layerCWTrange',np.arange(2,5))
    PBLCWTrange=kwargs.get('PBLCWTrange',np.arange(2,10))
    savemasks=kwargs.get('savemasks',True)
    savemaskname=kwargs.get('savemaskname','testmasksall.h5')
    
    mpltest = mtools.MPL()
    mpltest.fromHDF(filename)
    mpltest.time_resample(timestep=timestep)
    mpltest.calc_all()     
    if NRBmask:
        mpltest=NRB_mask_all(mpltest,NRBthreshold=NRBthresh)   
    
    molecular=molecular_detect(mpltest,varthresh=molthresh, winsize=molwinsize) 
    layers=find_layers(mpltest,noisethresh=layernoisethresh,cloudthresh=cloudthresh,
                       datatype=layerdtype,CWTwidth=CWTwidth,widths=layerCWTrange,
                       minwidth=layerminwidth)        
    mol_min=molecular.loc['Layer0']
    try:
        layer_min=layers.loc['Layer0']
    except KeyError:
        layer_min=mpltest.NRB[0].columns[-1]
    pbl=PBL_detect(mpltest,mol_min=mol_min,layer_min=layer_min,widths=PBLCWTrange)
    
    if savemasks:    
        store=pan.HDFStore(savemaskname)
        store['molecular']=molecular
        store['layers']=layers
        store['PBL']=pbl
        store.close()
    
    dictout={'mpl':mpltest,'molecular':molecular,'layers':layers,'pbl':pbl}
    return dictout   
    
def layermaskplot(mpl,**kwargs):
    fromHDF=kwargs.get('fromHDF',False)
    if fromHDF:
        HDFfile=kwargs.get('HDFfile',[])
        molecular=pan.read_hdf(HDFfile,'molecular')
        layers=pan.read_hdf(HDFfile,'layers')
        PBL=pan.read_hdf(HDFfile,'PBL')
    else:
        molecular=kwargs.get('molecular',[])
        layers=kwargs.get('layers',[])
        PBL=kwargs.get('PBL',[])
    SNRmask=kwargs.get('SNRmask',True)    
    SNRmasktype=kwargs.get('SNRmasktype','NRB')
    SNRthresh=kwargs.get('SNRthresh',1)
    hours=kwargs.get('hours',['03','06','09','12','15','18','21'])
    altrange=kwargs.get('altrange',np.arange(0,15030,30))
    saveplot=kwargs.get('saveplot',True)    
    plotfilepath=kwargs.get('plotfilepath',[])
    plotfilename=kwargs.get('plotfilename','testmaskplot.png')
    
    
    if SNRmask:
        MPLmasked = mpl.SNR[SNRmasktype][0]>=SNRthresh
        MPLmasked.replace(False,np.nan,inplace=True)
    
    mask,colordict=colormask_fromdict(mpl,PBL,molecular,layers)
    minalt=altrange[0]
    maxalt=altrange[-1]
    kwargs={'hours':hours,'altrange':(minalt,maxalt),'SNRmask':SNRmask,'saveplot':saveplot,
            'plotfilepath':plotfilepath,'plotfilename':plotfilename}
    colormask_plot(mask,colordict,**kwargs)

def scenemaker(layerdict,**kwargs):
    """
    Takes outputs from layer masking subroutines and generates a single pandas panel with all information
    needed to enter extinction processing phase
    
    Inputs:
    layerdict - dict objct containing the following key-value pairs
    
    kwargs:
    PBLrat - 
    
    """
    mpl=layerdict['mpl']
    molecular=layerdict['molecular']
    layers=layerdict['layers']
    PBL=layerdict['pbl']
    
    PBLrat=kwargs.get('PBLrat',30.0)
    molrat=kwargs.get('molrat',8.0*np.pi/3.0)
    moldepol=kwargs.get('moldepol',0.0035)
    udefrat=kwargs.get('udefrat',8.0*np.pi/3.0)
    udefdepol=kwargs.get('udefdepol',0.0035)
    savefile=kwargs.get('savefile',False)
    savefilename=kwargs.get('savefilename','test.h5')
    colordict=kwargs.get('colordict',[])
    
    if not colordict:
        colordict={'molecular':0,
                   'PBL':1,
                   'ice':2,
                   'water':3,
                   'mixed':4,
                   'dust':5,
                   'smoke':6,
                   'water_soluble':7,
                   'unidentified':8}
                   
    alts=mpl.data[0].columns
    times=mpl.data[0].index
    #initialize dataframes 
    typemask=pan.DataFrame(index=times,columns=alts)
    subtypemask=pan.DataFrame(index=times,columns=alts)
    lrat=pan.DataFrame(index=times,columns=alts)
    depolrat=pan.DataFrame(index=times,columns=alts)
    delta=pan.DataFrame(index=times,columns=alts)
    layerbase=pan.DataFrame(index=times,columns=alts)
    layertop=pan.DataFrame(index=times,columns=alts)
    colormask=pan.DataFrame(index=times,columns=alts)
    

    
    #fill dataframes with values defined by layerdict
    for t in times:
        tempmol=molecular.ix[:,t].dropna(axis=1,how='all')
        templayers=layers.ix[:,t].dropna(axis=1,how='all')
        tempPBL=PBL.ix[t]
        #first assign PBL props to altitudes below PBL height
        typemask.loc[t,typemask.columns<=tempPBL]='PBL'
        subtypemask.loc[t,subtypemask.columns<=tempPBL]='PBL'
        lrat.loc[t,lrat.columns<=tempPBL]=PBLrat
        depoltemp=mpl.depolrat[0].ix[t]
        depolrat.loc[t,depolrat.columns<=tempPBL]=np.mean(depoltemp[depoltemp.index<tempPBL])
        layerbase.loc[t,layerbase.columns<=tempPBL]=alts[0]
        layertop.loc[t,layertop.columns<=tempPBL]=tempPBL
        colormask.loc[t,colormask.columns<=tempPBL]=colordict['PBL']
        #then assign molecular props to altitudes identified as such
        for m in tempmol.iteritems():
            base=m[1]['Base']
            top=m[1]['Top']
            typemask.loc[t,(typemask.columns>=base)&(typemask.columns<=top)]='molecular'
            subtypemask.loc[t,(subtypemask.columns>=base)&(subtypemask.columns<=top)]='molecular'            
            lrat.loc[t,(lrat.columns>=base)&(lrat.columns<=top)]=molrat
            depolrat.loc[t,(depolrat.columns>=base)&(depolrat.columns<=top)]=moldepol
            delta.loc[t,(delta.columns>=base)&(delta.columns<=top)]=0.0
            layerbase.loc[t,(layerbase.columns>=base)&(layerbase.columns<=top)]=base
            layertop.loc[t,(layertop.columns>=base)&(layertop.columns<=top)]=top  
            colormask.loc[t,(colormask.columns>=base)&(colormask.columns<=top)]=colordict['molecular']
        #finally assign layer properties to each layer one by one
        for l in templayers.iteritems():
            base=l[1]['Base']
            top=l[1]['Top']
            temprat=l[1]['Lidar_Ratio']
            tempdelta=l[1]['Delta']
            temptype=l[1]['Type']
            tempsubtype=l[1]['Sub-Type']
            tempdepol=l[1]['Depol']
            typemask.loc[t,(typemask.columns>=base)&(typemask.columns<=top)]=temptype
            subtypemask.loc[t,(subtypemask.columns>=base)&(subtypemask.columns<=top)]=tempsubtype
            lrat.loc[t,(lrat.columns>=base)&(lrat.columns<=top)]=temprat
            depolrat.loc[t,(depolrat.columns>=base)&(depolrat.columns<=top)]=tempdepol
            delta.loc[t,(delta.columns>=base)&(delta.columns<=top)]=tempdelta
            layerbase.loc[t,(layerbase.columns>=base)&(layerbase.columns<=top)]=base
            layertop.loc[t,(layertop.columns>=base)&(layertop.columns<=top)]=top
            colormask.loc[t,(colormask.columns>=base)&(colormask.columns<=top)]=colordict[tempsubtype]
    
        #find and classify unidentified areas
        tempprof=typemask.loc[t]
        tempprof.fillna('unidentified',inplace=True)
        tempmask=pan.Series(alts,index=[v=='unidentified' for v in tempprof])
        tempgroups=tempmask.groupby(level=0) 
        #assuming uniform altitude steps
        altstep=alts[1]-alts[0]
        for g in tempgroups:
            if g[0]:
                tempalts=g[1]
                for key,alt in groupby(enumerate(tempalts), lambda (i,x):i-(x-tempalts.iloc[0])/altstep):                    
                    layeralt=map(operator.itemgetter(1),alt)
                    if len(layeralt)==1:
                        continue
                    else:
                        base=layeralt[0]
                        top=layeralt[-1]
                        typemask.loc[t,(typemask.columns>=base)&(typemask.columns<=top)]='unidentified'
                        subtypemask.loc[t,(subtypemask.columns>=base)&(subtypemask.columns<=top)]='unidentified'
                        lrat.loc[t,(lrat.columns>=base)&(lrat.columns<=top)]=udefrat
                        depolrat.loc[t,(depolrat.columns>=base)&(depolrat.columns<=top)]=udefdepol
                        delta.loc[t,(delta.columns>=base)&(delta.columns<=top)]=0
                        layerbase.loc[t,(layerbase.columns>=base)&(layerbase.columns<=top)]=base
                        layertop.loc[t,(layertop.columns>=base)&(layertop.columns<=top)]=top
                        colormask.loc[t,(colormask.columns>=base)&(colormask.columns<=top)]=colordict['unidentified']
                        
    typemask.fillna(method='ffill',inplace=True)
    subtypemask.fillna(method='ffill',inplace=True) 
    lrat.fillna(method='ffill',inplace=True)
    depolrat.fillna(method='ffill',inplace=True)
    layerbase.fillna(method='ffill',inplace=True)
    layertop.fillna(method='ffill',inplace=True)
    colormask.fillna(method='ffill',inplace=True)
    paneldict={'Type':typemask,'Sub-Type':subtypemask,'Lidar_Ratio':lrat,'Depol':depolrat,
               'Delta':delta,'Base':layerbase,'Top':layertop,'Color':colormask}
    panelout=pan.Panel.from_dict(paneldict)
    
    mpl.scenepanel=[panelout]

    if savefile:
        mpl.save_to_HDF(savefilename)        

    return mpl

def Tfromextinction(extinctprof): 
    alts=extinctprof.index
    tau=0.0
    oldz=alts[0]
    oldex=extinctprof.ix[oldz]
    for z in alts[1:]:
        deltaz=z-oldz
        tempex=(oldex+extinctprof.ix[z])/2.0
        tau+=-2*tempex*deltaz
        oldz=z
        oldex=extinctprof.ix[z]    
    Tout=np.exp(tau)
    return Tout
    
def Tfromprofile(prof,refprof=[],**kwargs):
    #requires profile containing at least [numvals] points on top and bottom
    #of molecular scattering free of aerosols
    wave=kwargs.get('wave',532.0)
    average=kwargs.get('average',True)
    numvals=kwargs.get('numvals',10)
    
    alts=prof.index.values
    if not refprof:
        refprof=itools.molprof(alts,wave=wave)['vals']*alts**2
    
    if average:
        profabove=np.mean(prof.iloc[-numvals:])
        profbelow=np.mean(prof.iloc[:numvals])
        refabove=np.mean(refprof.iloc[-numvals:])
        refbelow=np.mean(refprof.iloc[:numvals])
    else:
        profabove=prof.iloc[-1]
        profbelow=prof.iloc[0]
        refabove=prof.iloc[-1]
        refbelow=prof.iloc[0]
    
    normrat=refbelow/profbelow
    Tout=profabove*normrat/refabove
    return Tout

def Tmolecular(minalt,maxalt,refprof=[],wave=532.0):
    alts=np.arange(minalt,maxalt+1,1)
    if not refprof:
        refprof=itools.molprof(alts,wave=wave)['sigma_t']
    Tout=refprof.iloc[-1]/refprof.iloc[0]
    return Tout

def Tfrombounds(molbelow,molabove,wave=532.0):
    refalts=np.arange(molbelow.index[0],molabove.index[-1]+1,1)
    refprof=itools.molprof(refalts,wave=wave)['vals']*refalts**2
    refbelow=refprof.ix[molbelow.index]
    refabove=refprof.ix[molabove.index]        
    normrat=np.mean(refbelow.values)/np.mean(molbelow.values)
    Tout=np.mean(molabove.values)*normrat/np.mean(refabove.values)
    return Tout

def coefprofs(profin,lratprof,**kwargs):
    method=kwargs.get('method','klett2')
    refalt=kwargs.get('refalt',[])
    energy=kwargs.get('energy',1.0)
    calrange=kwargs.get('calrange',[])
    wave=kwargs.get('wave',532.0)
    
    if np.all(np.isnan(profin)):
        backprof=pan.Series(0,index=profin.index)
        extprof=backprof
    else:    
        if method=='fernald':
            if not calrange:
                calrange=[profin.index.values[-20],profin.index.values[-1]]
            lratin=np.mean(lratprof)
            backprof,extprof=itools.fernald(profin,lratin,E=energy,calrange=calrange)
        elif method=='klett':
            if not refalt:
                refalt=profin.index[-1]
            refext=itools.molprof([refalt],wave).loc['sigma_t'].values[0]
            backprof,extprof=itools.klett(profin,lratprof,r_m=refalt,sigma_m=refext)
        elif method=='klett2':
            lratin=np.mean(lratprof)
            backprof,extprof=itools.klett2(profin,lratprof)
    return backprof,extprof
    
def basiccorrection(mpl,**kwargs):    
    wave=kwargs.get('wave',532.0)
    refalt=kwargs.get('refalt',[])
    method=kwargs.get('method','klett2')
    lrat=kwargs.get('lrat',[])
    mode=kwargs.get('mode','copol')
    
    times=mpl.data[0].index
    alts=mpl.data[0].columns
    backscatter_copol=pan.DataFrame(data=np.nan,index=times,columns=alts,dtype='float')
    extinction_copol=pan.DataFrame(data=np.nan,index=times,columns=alts,dtype='float')
    backscatter_crosspol=pan.DataFrame(data=np.nan,index=times,columns=alts,dtype='float')
    extinction_crosspol=pan.DataFrame(data=np.nan,index=times,columns=alts,dtype='float')
    
    if mode=='copol' or mode=='both':
        NRB=mpl.NRB[0]
        scenepanel=mpl.scenepanel[0]
        
        for t in times:
            tempprof=NRB.ix[t]
            if not refalt:
                for a in alts[::-1]:
                    if not np.isnan(tempprof.ix[a]):
                        refalt=a
                        break
            if lrat:
                templrat=pan.Series(data=lrat,index=alts)
            else:
                templrat=scenepanel.loc['Lidar_Ratio',t]
            coefkwargs={}
            coefkwargs['wave']=wave
            coefkwargs['refalt']=refalt
            coefkwargs['method']=method
            coefkwargs['energy']=mpl.header['energy'].ix[t]
            coefkwargs['calrange']=[a-100,a]
            tempbeta,tempsigma=coefprofs(tempprof,templrat,**coefkwargs)
            backscatter_copol.ix[t]=tempbeta
            extinction_copol.ix[t]=tempsigma
    
    if mode=='crosspol' or mode=='both':
        NRB=mpl.NRB[1]
        scenepanel=mpl.scenepanel[1]        
        for t in times:
            tempprof=NRB.ix[t]
            if not refalt:
                for a in alts[::-1]:
                    if not np.isnan(tempprof.ix[a]):
                        refalt=a
                        break
            if lrat:
                templrat=pan.Series(data=lrat,index=alts)
            else:
                templrat=scenepanel.loc['Lidar_Ratio',t]
            coefkwargs={}
            coefkwargs['wave']=wave
            coefkwargs['refalt']=refalt
            coefkwargs['method']=method
            coefkwargs['energy']=mpl.header['energy'].ix[t]
            coefkwargs['calrange']=[a-100,a]
            tempbeta,tempsigma=coefprofs(tempprof,templrat,**coefkwargs)
            backscatter_crosspol.ix[t]=tempbeta
            extinction_crosspol.ix[t]=tempsigma
        
    mpl.backscatter=[backscatter_copol,backscatter_crosspol]        
    mpl.extinction=[extinction_copol,extinction_crosspol]
        
    return mpl

def Tmatcher(lrat,profin,Tin,method,layerrange,wave,energy,refalt,calrange):    
    alts=profin.index
    if not refalt:
        for a in alts[::-1]:
            if not np.isnan(profin.ix[a]):
                refalt=a
                break
    
    coefkwargs={}
    coefkwargs['wave']=wave
    coefkwargs['refalt']=refalt
    coefkwargs['method']=method
    coefkwargs['calrange']=calrange
    coefkwargs['energy']=energy
    lratprof=pan.Series(data=lrat,index=profin.index)
    backprof,extprof=coefprofs(profin,lratprof,**coefkwargs)    
    Tcalc=Tfromextinction(extprof.ix[layerrange[0]:layerrange[-1]])
    return Tin-Tcalc

def lratpart(beta_t,lrat_t,wave=532.0,tol=1e-2):
    lrat_m=8.0*np.pi/3.0
    
    alts=beta_t.index.values
    beta_m=itools.molprof(alts,wave=wave).loc['beta_t']
    
    diffbeta=(beta_t-beta_m)/beta_m
    lrat_p=(lrat_t*beta_t-lrat_m*beta_m)/(beta_t-beta_m)
    lrat_out=lrat_p*diffbeta.gt(tol)
    lrat_out.fillna(0.0,inplace=True)
    
    return lrat_out

def simplelayercorrection(layerprof,backprof,extprof,lrat0,**kwargs):
    #profin must be bounded 
    molbelow=kwargs.get('molbelow',[])
    molabove=kwargs.get('molabove',[])
    method=kwargs.get('method','klett2')
    lratrange=kwargs.get('lratrange',[10,100])
    tolerance=kwargs.get('tolerance',0.01)
    numreps=kwargs.get('numreps',100)
    wave=kwargs.get('wave',532.0)
    average=kwargs.get('average',True)
    numvals=kwargs.get('numvals',10)
    partrat=kwargs.get('partrat',True)
    energy=kwargs.get('energy',1.0)
    #calculate optical depth from profile or boundaries
    if len(molbelow)>0 and len(molabove)>0:
        Tin=Tfrombounds(molbelow,molabove)
    else:
        Tin=Tfromprofile(layerprof,wave=wave,average=average,numvals=numvals)    
    #calculate optical depth from lidar ratio
    layerrange=[layerprof.index[0],layerprof.index[-1]]    
    if method=='klett2':
        refalt=[]
        calrange=[]
        try:
            lratcor=opt.newton(Tmatcher,lrat0,args=(layerprof,Tin,method,layerrange,wave,energy,refalt,calrange),maxiter=numreps,tol=tolerance)
            print "{0} successfully corrected to {1} after < {2} iterations".format(lrat0,lratcor,numreps)
        except RuntimeError:
            lratcor=lrat0
            print "Warning! Failed to converge after {0} iterations".format(numreps)
        coefkwargs={}
        coefkwargs['refalt']=refalt
        coefkwargs['method']=method
        coefkwargs['wave']=wave
        coefkwargs['calrange']=calrange
        coefkwargs['energy']=energy
        if lratrange[0]<=lratcor<=lratrange[1]:
            lratprof=pan.Series(data=lratcor,index=layerprof.index)
            backprof,extprof=coefprofs(layerprof,lratprof,**coefkwargs)        
        else:
            lratprof=pan.Series(data=lrat0,index=layerprof.index)
            backprof,extprof=coefprofs(layerprof,lratprof,**coefkwargs) 
        if partrat:
                lratout=lratpart(backprof,lratcor)
                lratprof=pan.Series(data=lratout,index=layerprof.index) 
    else:
        tempprof=layerprof.append(molabove)
        refalt=tempprof.index[-1]
        if len(molabove)>numvals:
            calrange=[tempprof.index[-numvals],tempprof.index[-1]]
        else:
            calrange=[molabove.index[0],molabove.index[-1]]
        try:
            lratcor=opt.newton(Tmatcher,lrat0,args=(tempprof,Tin,method,layerrange,wave,energy,refalt,calrange),maxiter=numreps,tol=tolerance)
            print "{0} successfully corrected to {1} after < {2} iterations".format(lrat0,lratcor,numreps)
        except RuntimeError:
            lratcor=lrat0
            print "Warning! Failed to converge after {0} iterations".format(numreps)
        coefkwargs={}
        coefkwargs['refalt']=refalt
        coefkwargs['method']=method
        coefkwargs['wave']=wave
        coefkwargs['calrange']=calrange
        coefkwargs['energy']=energy
        if lratrange[0]<=lratcor<=lratrange[1]:
            lratout=lratcor
        else:
            lratout=lrat0
        molrat=8.0*np.pi/3.0
        lratprof_ref=pan.Series(data=molrat,index=molabove.index)
        lratprof=pan.Series(data=lratout,index=layerprof.index)
        templratprof=lratprof.append(lratprof_ref)
        tempbackprof,tempextprof=coefprofs(tempprof,templratprof,**coefkwargs)
        backprof=tempbackprof.ix[layerprof.index[0]:layerprof.index[-1]]
        extprof=tempextprof.ix[layerprof.index[0]:layerprof.index[-1]]
        if partrat:
            lratprof=lratpart(backprof,lratout)
    
    return backprof,extprof,lratprof

def sceneoptimize(mplin,**kwargs):
    wave=kwargs.get('wave',532.0)
    refalt=kwargs.get('refalt',[])
    method=kwargs.get('method','klett2')
    calrange=kwargs.get('calrange',[])
    lrat=kwargs.get('lrat',[])
    inplace=kwargs.get('inplace',True)
    mode=kwargs.get('mode','copol')
        
    NRB=mplin.NRB[0]
    times=NRB.index
    alts=NRB.columns
    if inplace:
        mplout=mplin
    else:
        mplout=deepcopy(mplin)
    
    #start by performing initial backscatter and extinction corrections using initial lrat values
    coefkwargs={}
    coefkwargs['refalt']=refalt
    coefkwargs['method']=method
    coefkwargs['wave']=wave
    coefkwargs['calrange']=calrange
    coefkwargs['lrat']=lrat
    coefkwargs['mode']=mode
    
    if len(mplout.backscatter)==0 or len(mplout.extinction)==0:
        mplout=basiccorrection(mplout,**coefkwargs)
        
    if mode=='copol' or mode=='both':
        backscatter=mplout.backscatter[0]
        extinction=mplout.extinction[0]
        try:
            panelout=mplout.scenepanel[0]
        except IndexError:
            print "Missing {0} Scene Panel".format(mode)
            return mplout
        #then go through profiles one by one and perform corrections for each layer
        dflrat=panelout.loc['Lidar_Ratio']
        for t in times:
            NRBprof=NRB.loc[t]
            lratprof=panelout.loc['Lidar_Ratio',t]
            baseprof=panelout.loc['Base',t]
            topprof=panelout.loc['Top',t]
            typeprof=panelout.loc['Type',t]
            backprof=backscatter.loc[t]
            extprof=extinction.loc[t]
            a=alts[0]
            while True:
                layerbase=baseprof[a]
                layertop=topprof[a]
                lrat=lratprof[a]
                ltype=typeprof[a]
                lratbelow=lratprof[lratprof.index<layerbase]
                lratabove=lratprof[lratprof.index>layertop]
                typebelow=typeprof[typeprof.index<layerbase]
                typeabove=typeprof[typeprof.index>layertop]
                if len(lratabove)==0:
                    break
                if len(lratbelow)==0:
                    a=lratabove.index[0]
                    continue
                if ltype=='PBL' or ltype=='molecular' or ltype=='unidentified':
                    a=lratabove.index[0]
                else:                
                    if typebelow.iloc[-1]=='molecular' and typeabove.iloc[0]=='molecular':
                        tempNRB=NRBprof[(NRBprof.index>=layerbase)&(NRBprof.index<=layertop)]
                        tempback=backprof[(backprof.index>=layerbase)&(backprof.index<=layertop)]
                        tempext=extprof[(extprof.index>=layerbase)&(extprof.index<=layertop)]
                        tempprof=[]
                        z=[]
                        for tempalt in lratbelow.index[::-1]:
                            if typebelow.ix[tempalt]=='molecular':
                                z.append(tempalt)
                                tempprof.append(NRBprof.ix[tempalt])
                            else:
                                break
                        molbelow=pan.Series(tempprof[::-1],index=z[::-1])
                        tempprof=[]
                        z=[]
                        for tempalt in lratabove.index:
                            if typeabove.ix[tempalt]=='molecular':
                                z.append(tempalt)
                                tempprof.append(NRBprof.ix[tempalt])
                            else:
                                break
                        molabove=pan.Series(tempprof,index=z)
                        corkwargs={}
                        corkwargs['energy']=mplin.header['energy'].ix[t]
                        corkwargs['method']=method
                        corkwargs['molbelow']=molbelow
                        corkwargs['molabove']=molabove
                        
                        print "Type is {0} mean signal is {1}".format(ltype,np.mean(tempNRB))
                        backcor,extcor,lratcor=simplelayercorrection(tempNRB,tempback,tempext,lrat,**corkwargs)
                        backprof[(backprof.index>=layerbase)&(backprof.index<=layertop)]=backcor
                        extprof[(extprof.index>=layerbase)&(extprof.index<=layertop)]=extcor
                        lratprof[(lratprof.index>=layerbase)&(lratprof.index<=layertop)]=lratcor
                
                    a=lratabove.index[0]
            backscatter.loc[t]=backprof
            extinction.loc[t]=extprof
            dflrat.loc[t]=lratprof
        panelout.loc['Lidar_Ratio']=dflrat
        
        mplout.backscatter[0]=backscatter
        mplout.extinction[0]=extinction
        mplout.scenepanel[0]=panelout
    
    if mode=='crosspol' or mode=='both':
        backscatter=mplout.backscatter[1]
        extinction=mplout.extinction[1]
        try:
            panelout=mplout.scenepanel[1]
        except IndexError:
            print "Missing {0} Scene Panel".format(mode)
            return mplout
        #then go through profiles one by one and perform corrections for each layer
        dflrat=panelout.loc['Lidar_Ratio']
        for t in times:
            NRBprof=NRB.loc[t]
            lratprof=panelout.loc['Lidar_Ratio',t]
            baseprof=panelout.loc['Base',t]
            topprof=panelout.loc['Top',t]
            typeprof=panelout.loc['Type',t]
            backprof=backscatter.loc[t]
            extprof=extinction.loc[t]
            a=alts[0]
            while True:
                layerbase=baseprof[a]
                layertop=topprof[a]
                lrat=lratprof[a]
                ltype=typeprof[a]
                lratbelow=lratprof[lratprof.index<layerbase]
                lratabove=lratprof[lratprof.index>layertop]
                typebelow=typeprof[typeprof.index<layerbase]
                typeabove=typeprof[typeprof.index>layertop]
                if len(lratabove)==0:
                    break
                if len(lratbelow)==0:
                    a=lratabove.index[0]
                    continue
                if ltype=='PBL' or ltype=='molecular' or ltype=='unidentified':
                    a=lratabove.index[0]
                else:                
                    if typebelow.iloc[-1]=='molecular' and typeabove.iloc[0]=='molecular':
                        tempNRB=NRBprof[(NRBprof.index>=layerbase)&(NRBprof.index<=layertop)]
                        tempback=backprof[(backprof.index>=layerbase)&(backprof.index<=layertop)]
                        tempext=extprof[(extprof.index>=layerbase)&(extprof.index<=layertop)]
                        tempprof=[]
                        z=[]
                        for tempalt in lratbelow.index[::-1]:
                            if typebelow.ix[tempalt]=='molecular':
                                z.append(tempalt)
                                tempprof.append(NRBprof.ix[tempalt])
                            else:
                                break
                        molbelow=pan.Series(tempprof[::-1],index=z[::-1])
                        tempprof=[]
                        z=[]
                        for tempalt in lratabove.index:
                            if typeabove.ix[tempalt]=='molecular':
                                z.append(tempalt)
                                tempprof.append(NRBprof.ix[tempalt])
                            else:
                                break
                        molabove=pan.Series(tempprof,index=z)
                        corkwargs={}
                        corkwargs['energy']=mplin.header['energy'].ix[t]
                        corkwargs['method']=method
                        corkwargs['molbelow']=molbelow
                        corkwargs['molabove']=molabove
                        
                        print "Type is {0} mean signal is {1}".format(ltype,np.mean(tempNRB))
                        backcor,extcor,lratcor=simplelayercorrection(tempNRB,tempback,tempext,lrat,**corkwargs)
                        backprof[(backprof.index>=layerbase)&(backprof.index<=layertop)]=backcor
                        extprof[(extprof.index>=layerbase)&(extprof.index<=layertop)]=extcor
                        lratprof[(lratprof.index>=layerbase)&(lratprof.index<=layertop)]=lratcor
                
                    a=lratabove.index[0]
            backscatter.loc[t]=backprof
            extinction.loc[t]=extprof
            dflrat.loc[t]=lratprof
        panelout.loc['Lidar_Ratio']=dflrat
        
        mplout.backscatter[1]=backscatter
        mplout.extinction[1]=extinction
        mplout.scenepanel[1]=panelout
        
    return mplout

def quickplot(df,**kwargs):
    ar=kwargs.get('ar',2.0)
    vrange=kwargs.get('vrange',[0,1])
    fsize=kwargs.get('fsize',21)
    maptype=kwargs.get('maptype','customjet')
    orientation=kwargs.get('orientation','Vertical')
    overcolor=kwargs.get('overcolor','w')
    undercolor=kwargs.get('undercolor','k')
    numvals=kwargs.get('numvals',50)
    
    fig=plt.figure()
    ax=fig.add_subplot(111)    
    data=df.values.T[::-1]
    xdata=df.index
    ydata=df.columns[::-1]
    
    my_cmap=mplot.custom_cmap(maptype=maptype,numvals=numvals,overcolor=overcolor,
                              undercolor=undercolor)  
    im = ax.imshow(data,vmin=vrange[0],vmax=vrange[1],cmap = my_cmap)
    mplot.forceAspect(ax,ar)       
    mplot.altticks(ax, ydata, fsize = fsize, tcolor = 'w')
    
    if orientation=='vertical':
        ax.set_ylabel('Altitude [m]', fontsize = fsize+4, labelpad = 15)
    elif orientation=='vorizontal':
        ax.set_ylabel('Horizontal Range [m]', fontsize = fsize+4, labelpad = 15)
    for line in ax.yaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(1)        
    ax.axis('tight')
    vstep=(vrange[1]-vrange[0])/5.0
    cbar = fig.colorbar(im, orientation = 'vertical', aspect = 6, extend='both')
    cbar.set_ticks(np.arange(vrange[0],vrange[1]+vstep,vstep))
    cbar.set_ticklabels(np.arange(vrange[0],vrange[1]+vstep,vstep))
    cbar.ax.tick_params(labelsize = fsize)
    cbar.ax.set_ylabel('$[counts*km^{2}/(\mu s*\mu J)$')
    mplot.dateticks(ax, xdata, hours = ['03','06', '09','12', '15','18','21'], 
              fsize = fsize, tcolor = 'w')
    fig.canvas.draw()
    
#def multifile_sceneproc():
#        os.chdir('C:\Users\dashamstyr\Dropbox\Lidar Files\UBC Cross-Cal\\20131014-20131016\\10-15\Processed')
#        plotfilepath='C:\Users\dashamstyr\Dropbox\Lidar Files\UBC Cross-Cal\\20131014-20131016\\10-15\Figures'
#        filepath=mtools.get_files('Select processed file',filetype=('.h5','*.h5'))[0]
#        filename=os.path.split(filepath)[-1]
#        plotfilename='{0}_coefplot.png'.format(filename.split('_proc')[0])
#        timestep='120S'
#        SNRthreshold=3
#        layerdict=findalllayers(filename,timestep=timestep,NRBmask=False)
#        scenepanel=scenemaker(layerdict)
#        mpltemp=layerdict['mpl']
#        backscatter,extinction=basiccorrection(mpltemp,scenepanel)
#        mpltemp.backscatter.append(backscatter)
#        mpltemp.extinction.append(extinction)
#        
#        altrange = np.arange(150,15030,30)
#        starttime = []
#        endtime = []
#        timestep = '120S'
#        hours = ['03','06','09','12','15','18','21']
#        topplot_limits=(0.0,1e-7,2e-8)
#        bottomplot_limits=(0.0,2e-7,4e-8)
#        
#        kwargs = {'saveplot':True,'showplot':True,'verbose':True,
#                    'savefilepath':plotfilepath,'savefilename':plotfilename,
#                    'hours':hours,'bottomplot_limits':bottomplot_limits,'timestep':timestep,
#                    'topplot_limits':topplot_limits,'SNRmask':False,'altrange':altrange,
#                    'colormap':'customjet','orientation':'vertical','toptype':'backscatter',
#                    'bottomtype':'extinction'}
#        
#        mplot.doubleplot(mpltemp,**kwargs)    
#        
#        layermaskplot(mpl=layerdict['mpl'],molecular=layerdict['molecular'],layers=layerdict['layers'],
#                      PBL=layerdict['pbl'],saveplot=True,plotfilepath=plotfilepath,plotfilename=plotfilename)

    
if __name__=='__main__':
    os.chdir('C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Aksu Data')
    savefilepath='C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Aksu Data\Processed'
    plotfilepath='C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Aksu Data\Figures'
    filepath=mtools.get_files('Select processed file',filetype=('.h5','*.h5'))[0]
    filename=os.path.split(filepath)[-1]
    savefilename='{0}_scenepanel.h5'.format(filename.split('_proc')[0])
    plotfilename='{0}_coefplot.png'.format(filename.split('_proc')[0])
    timestep='120S'
    SNRthreshold=0.1
    molthresh=0.1
    layernoisethresh=0.1
    layerdict=findalllayers(filepath,timestep=timestep,NRBmask=False,molthresh=molthresh,
                            layernoisethresh=layernoisethresh)
    mpltemp=scenemaker(layerdict)
    mpltemp=basiccorrection(mpltemp,inplace=True)
#    mpltemp.save_to_HDF(savefilename)
#    mplemp2=sceneoptimize(mpltemp,method='klett2',inplace=False)
    
    altrange = np.arange(150,9030,30)
    starttime = []
    endtime = []
    timestep = '120S'
    hours = ['03','06','09','12','15','18','21']
    topplot_limits=(0.0,1e-7,2e-8)
    bottomplot_limits=(0.0,2e-6,4e-7)
    
    kwargs = {'saveplot':True,'showplot':True,'verbose':True,
                'savefilepath':plotfilepath,'savefilename':plotfilename,
                'hours':hours,'bottomplot_limits':bottomplot_limits,'timestep':timestep,
                'topplot_limits':topplot_limits,'SNRmask':False,'altrange':altrange,
                'colormap':'customjet','orientation':'vertical','toptype':'backscatter',
                'bottomtype':'extinction'}
    
    mplot.doubleplot(mpltemp,**kwargs)    
    
    layermaskplot(mpl=layerdict['mpl'],molecular=layerdict['molecular'],layers=layerdict['layers'],
                  PBL=layerdict['pbl'],saveplot=True,plotfilepath=plotfilepath,plotfilename=plotfilename)
#    mpl=mtools.MPL()
#    mpl.fromHDF(filename)
#    mpl.time_resample(timestep=timestep)
#    mpl.calc_all()
    
    #layermaskplot(mpltemp,fromHDF=True,HDFfile='testmasksall.h5')
    
#    mpltest=mtools.MPL()
#    mpltest.fromHDF(filename)
#    mpltest.time_resample(timestep=timestep)
#    mpltest.calc_all()
#    
#    mplmasked=SNR_mask_all(mpltest,SNRthreshold=SNRthreshold)
#    
#    
##    maxalt=9000    
##    NRBprof=mpltest.NRB[0].ix['2014-05-03T00:00:00.000000000-0700']
##    NRBprof=NRBprof[NRBprof.index<maxalt]
##    dataprof=mpltest.data[0].ix['2014-05-03T00:00:00.000000000-0700']
##    dataprof=dataprof[dataprof.index<maxalt]
##    molecular_detect_single(NRBprof,dataprof)
#    
#    molecular=molecular_detect(mpltest,varthresh=1, winsize=5) 
#    layers=find_layers(mpltest,noisethresh=0.4,cloudthresh=(1.0,0.20),datatype='NRB',layerwidth=4,widths=np.arange(2,5))
#    
##    with pan.get_store('testmolecular.h5') as molstore:
##        molecular=molstore['molecular']
##    
##    with pan.get_store('testlayers.h5') as laystore:
##        layers=laystore['layers'] 
#        
#    mol_min=molecular.loc['Layer0']
#    layer_min=layers.loc['Layer0']
#    
##    profile=mpltest.NRB[0].iloc[20,:]
##    l=layers.iloc[:,20,:]
##    
##    layerprofplot(profile,l,numlayer=10)
#    pbl=PBL_detect(mpltest,mol_min=mol_min,layer_min=layer_min,widths=np.arange(2,10))
#    
#    
#    SNRmask = mpltest.SNR['data'][0]>=SNRthreshold
#    SNRmask.replace(False,np.nan,inplace=True)
#    
#    mask,colordict=colormask_fromdict(mpltest,pbl,molecular,layers)
#    
#    hours=['03','06','09','12','15','18','21']
#    altrange=altrange = np.arange(0,15030,30)
#    minalt=altrange[0]
#    maxalt=altrange[-1]
#    colormask_plot(mask,colordict,hours=hours,altrange=(minalt,maxalt),SNRmask=SNRmask)
#    
##    mpltest.alt_resample(altrange)
##    mpltest.calc_all() 
##    kwargs = {'saveplot':False,'showplot':True,'verbose':True,'altrange':altrange}    
##    mplot.doubleplot(mplmasked,**kwargs)
#
#    proftime='2014-05-02T23:30:00.000000000-0700'
#    
#    prof=mpltest.NRB[0].ix[proftime]
#    layersin=layers.ix[:,proftime,:]
#    layerprofplot(prof,layersin,numlayer=30)