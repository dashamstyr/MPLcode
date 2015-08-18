# -*- coding: utf-8 -*-
"""
Created on Tue May 19 12:50:39 2015

@author: User
"""
import os,sys,glob

import numpy as np
import MPLtools as mtools
import MPLprocesstools as mproc
import MPLfileproc as mfile
import pandas as pan
import datetime
import pytz
import ephem
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns

def PBLanalyze(**kwargs):
    """
    Program used to find and hsitogram PBL heights and capping layers
    """
    filedir=kwargs.get('filedir',[])
    location=kwargs.get('location',(48.9167,-125.5333))
    altrange = kwargs.get('altrange',np.arange(0.150,15.030,0.030))        
    timestep = kwargs.get('timestep','60S')
    molthresh=kwargs.get('molthresh',1.0)
    layernoisethresh=kwargs.get('layernoisethresh',1.0)      
    bg_alt=kwargs.get('bg_alt',None)
    datatype=kwargs.get('datatype','NRB')
    winsize=kwargs.get('winsize',5)
    wavelet=kwargs.get('wavelet',signal.ricker)
    noisethresh=kwargs.get('noisethresh',0.4)
    cloudthresh=kwargs.get('cloudthresh',(1.0,0.20))
    CWTwidth=kwargs.get('CWTwidth',2)
    minwidth=kwargs.get('minwidth',4)
    layerCWTrange=kwargs.get('layerCWTrange',np.arange(2,5))
    PBLwavelet=kwargs.get('PBLwavelet',mproc.dog)
    PBLCWTrange=kwargs.get('PBLCWTrange',np.arange(2,10))
    sigma0=kwargs.get('sigma0',0.01)
    waterthresh=kwargs.get('waterthresh',0.10)
    icethresh=kwargs.get('icethresh',0.35)
    smokethresh=kwargs.get('smokethresh',0.10)
    dustthresh=kwargs.get('dustthresh',0.20)
    
    if not filedir:
        filedir=mtools.set_dir('Select file directory to process')
    
    olddir=os.getcwd()
    os.chdir(filedir)    
    rawfiles = glob.glob('*.mpl')
    rawfiles.sort()
    
    [path,startfile] = os.path.split(rawfiles[0])
    [path,endfile] = os.path.split(rawfiles[-1])
    startdate = kwargs.get('startdate',mfile.MPLtodatetime(startfile))
    enddate = kwargs.get('enddate',mfile.MPLtodatetime(endfile))
    mpllist=[]    
    datelist=[]
    PBL=pan.DataFrame()
    caplayer=pan.DataFrame()
    layerkwargs={'altrange':altrange,'timestep':timestep,'bg_alt':bg_alt,'datatype':datatype,
             'molthresh':molthresh,'winsize':winsize,'layernoisethresh':layernoisethresh,
             'wavelet':wavelet,'noisethresh':noisethresh,'cloudthresh':cloudthresh,
             'CWTwidth':CWTwidth,'minwidth':minwidth,'layerCWTrange':layerCWTrange,
             'PBLwavelet':PBLwavelet,'PBLCWTrange':PBLCWTrange,'sigma0':sigma0,
             'waterthresh':waterthresh,'icethresh':icethresh,'smokethresh':smokethresh,
             'dustthresh':dustthresh}
    
               
    filedates=[mfile.MPLtodatetime(os.path.split(f)[1]) for f in rawfiles]
    filteredfiles=[(r,t) for r,t in zip(rawfiles,filedates) if (t>=startdate) and (t<=enddate)] 
    
    oldday=filedates[0]
    tempephem,tempperiods=get_ephemera(location,oldday)
    tempmpl=[]
    tempdate=[]
    ephemera=[]
    periods=[]
    for f,d in filteredfiles:
        newday=d
        if newday.date()==oldday.date():
            tempmpl.append(f)
            tempdate.append(d)
        else:
            mpllist.append(tempmpl)
            datelist.append(tempdate)
            ephemera.append(tempephem)
            periods.append(tempperiods)
            tempmpl=[f]
            tempdate=[d]
            tempephem,tempperiods=get_ephemera(location,newday)
        if f==filteredfiles[-1][0]:
            mpllist.append(tempmpl)
            datelist.append(tempdate)
            ephemera.append(tempephem)
            periods.append(tempperiods)
        oldday=newday
            
    
    for mplgrp,datgrp,e,p in zip(mpllist,datelist,ephemera,periods):    
        for n in range(len(e)-1):
            tempmpl=[m for m,t in zip(mplgrp,datgrp) if (t>=e[n])and(t<e[n+1])]
            if len(tempmpl)==0:
                continue
            else:                
                templayerdict=get_layers(tempmpl,**layerkwargs)
                tempPBL=pan.DataFrame(data=templayerdict['pbl'],columns=[p[n]])
                tempcap=find_capper(templayerdict,p[n])
                PBL=PBL.append(tempPBL)
                caplayer=caplayer.append(tempcap)
                
    os.chdir(olddir)
    return PBL,caplayer
    

def get_ephemera(loc,datein,localzone=pytz.timezone('US/Pacific')):
    daystart=datetime.datetime(datein.year,datein.month,datein.day,0,0)
    daystart_local=localzone.localize(daystart)
    daystart_utc=daystart_local.astimezone(pytz.utc)
    dayend=daystart+datetime.timedelta(days=1)
    o=ephem.Observer()
    o.lat=str(loc[0])
    o.long=str(loc[1])
    sun=ephem.Sun()
    sunrise=ephem.localtime(o.next_rising(sun,start=daystart_utc))
    solarnoon=ephem.localtime(o.next_transit(sun,start=daystart_utc))
    sunset=ephem.localtime(o.next_setting(sun,start=daystart_utc))
    solarmidnight=ephem.localtime(o.next_antitransit(sun,start=daystart_utc))
    
    if solarmidnight>sunrise:
        periods=['Night','Morning','Afternoon','Evening','Night']        
        ephemera=[daystart,sunrise,solarnoon,sunset,solarmidnight,dayend]
    else:
        periods=['Evening','Night','Morning','Afternoon','Evening']        
        ephemera=[daystart,solarmidnight,sunrise,solarnoon,sunset,dayend]
    return ephemera,periods

def get_layers(rawfiles,**kwargs):
    altrange = kwargs.get('altrange',np.arange(0.150,15.03,0.030))        
    timestep = kwargs.get('timestep','60S')
    molthresh=kwargs.get('molthresh',1.0)
    layernoisethresh=kwargs.get('layernoisethresh',1.0)      
    bg_alt=kwargs.get('bg_alt',None)
    datatype=kwargs.get('datatype','NRB')
    winsize=kwargs.get('winsize',5)
    wavelet=kwargs.get('wavelet',signal.ricker)
    noisethresh=kwargs.get('noisethresh',0.4)
    cloudthresh=kwargs.get('cloudthresh',(1.0,0.20))
    CWTwidth=kwargs.get('CWTwidth',2)
    minwidth=kwargs.get('minwidth',4)
    layerCWTrange=kwargs.get('layerCWTrange',np.arange(2,5))
    PBLwavelet=kwargs.get('PBLwavelet',mproc.dog)
    PBLCWTrange=kwargs.get('PBLCWTrange',np.arange(2,10))
    sigma0=kwargs.get('sigma0',0.01)
    waterthresh=kwargs.get('waterthresh',0.10)
    icethresh=kwargs.get('icethresh',0.35)
    smokethresh=kwargs.get('smokethresh',0.10)
    dustthresh=kwargs.get('dustthresh',0.20)
    starttime=kwargs.get('starttime',datetime.datetime(1970,1,1))
    endtime=kwargs.get('endtime',datetime.datetime.today())    
    
#    MPLdat_event=mtools.MPL()
    for r in rawfiles:
        [path,tempname] = os.path.split(r)
        if starttime <= mfile.MPLtodatetime(tempname) <= endtime: 
            MPLdat_temp = mtools.MPL()
            MPLdat_temp.fromMPL(tempname)
            MPLdat_temp.alt_resample(altrange)    
            try:
                MPLdat_event.append(MPLdat_temp)
            except NameError:
                MPLdat_event = MPLdat_temp
       
    #sort by index to make certain data is in order then set date ranges to match
    MPLdat_event.header.sort_index()
    
    for n in range(MPLdat_event.header['numchans'][0]):
        data = MPLdat_event.data[n]
        data.sort_index()
    
    MPLdat_event.time_resample(timestep)
    MPLdat_event.calc_all()
    
    layerkwargs={'timestep':timestep,'bg_alt':bg_alt,'datatype':datatype,
                 'molthresh':molthresh,'winsize':winsize,'layernoisethresh':layernoisethresh,
                 'wavelet':wavelet,'noisethresh':noisethresh,'cloudthresh':cloudthresh,
                 'CWTwidth':CWTwidth,'minwidth':minwidth,'layerCWTrange':layerCWTrange,
                 'PBLwavelet':PBLwavelet,'PBLCWTrange':PBLCWTrange,'sigma0':sigma0,
                 'waterthresh':waterthresh,'icethresh':icethresh,'smokethresh':smokethresh,
                 'dustthresh':dustthresh}
                 
    layerdict=mproc.findalllayers(mplin=MPLdat_event,**layerkwargs)
    
    return layerdict

def find_capper(layerdict,dayperiod,typemode='Sub-Type'):
    molalt=layerdict['molecular']['Layer0']['Base']
    layeralt=layerdict['layers']['Layer0']['Base']
    layertype=layerdict['layers']['Layer0'][typemode]
    
    captype=pan.DataFrame(index=molalt.index,columns=[dayperiod])
    for i in captype.index:
        if layeralt.ix[i]<molalt.ix[i]:
            captype.ix[i]=layertype.ix[i]
        else:
            captype.ix[i]='Clear Air'
    
    return captype
        
#def pieplot(dfin,**kwargs):
#    
#    numfigs=len(plt.get_fignums())
#    if type(dfin.values[0,0])==str:
#        dfplot=dfin.apply(pan.value_counts)
#    else:
#        dfplot=dfin
#    
#    labels=dfplot.index.values
#    fig=plt.figure(numfigs+1)
#    subplot_counter=1
#    for c in dfplot.columns:
#        subplot_label=c
#        try:
#            ax=fig.add_subplot(layout[0],layout[1],subplot_counter)
#        except NameError:
#            ax=fig.add_subplot(1,len(dfplot.columns)+1,subplot_counter)
#
#        ax.pie(dfplot[c].values,autopct='%1.1f%%',radius=2,pctdistance=2.2)
#        ax.set_title(c)
#        ax.axis('equal')
#        subplot_counter+=1
#    
#    return fig
        
def barplot(dfin,**kwargs):
    scaled=kwargs.get('scaled',True)
    stacked=kwargs.get('stacked',True)
    numfigs=len(plt.get_fignums())
    
    if scaled==True:
        dfplot=100.0*dfin.div(dfin.sum(axis=1), axis=0)
        fmt='%.0f%%'
    else:
        dfplot=dfin
        fmt='%.00f%'
        
    sns.set_palette("Set3", len(dfplot.columns))
    sns.set_context('poster')
    sns.set_style('darkgrid',{'legend.frameon':True})
    ax=dfplot.plot(kind='barh',stacked=stacked)
    
    xticks=mtick.FormatStrFormatter(fmt)
    ax.xaxis.set_major_formatter(xticks)
    
def violinplot(dfin,**kwargs):
    
    sns.violinplot(dfin)
#    dfplot = [v.dropna() for k, v in dfin.iteritems()]
#    ax=sns.violinplot(dfplot) 
#    
#    xlabels=dfin.columns.values
#    xpos=range(len(xlabels))
#    plt.xticks(xpos,xlabels)
#    ax.set_axis_labels(xlabels)


if __name__=='__main__':
    
    os.chdir('C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\UBC Data')
    timestep='240S'
    altrange=np.arange(0.150,2.0,0.03)
#    startdate=datetime.datetime(2013,8,27,0)
#    enddate=datetime.datetime(2015,1,6,23)
    PBL,caplayer=PBLanalyze(timestep=timestep,altrange=altrange)
    store=pan.HDFStore('.\Processed\UBC_Boundary_Layer_stats_all.h5')
    store['PBL']=PBL
    store['caplayer']=caplayer
#    PBL=store['PBL']
#    caplayer=store['caplayer']
    
    store.close()
     
    violinplot(PBL)
    capplot=caplayer.apply(pan.value_counts).transpose()
    barplot(capplot)
#    piefig=pieplot(caplayer)    
    