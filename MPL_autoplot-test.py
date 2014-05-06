# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 17:04:59 2014

@author: dashamstyr

Program to generate processed image of files from most recent 24 hour window in miniMPL data file.

"""

import os, sys
import pandas as pan
import datetime,time
import numpy as np
import glob
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


codelib = '/data/lv1/pcottle/MPLCode'
datalib = '/data/lv1/pcottle/MPLData/Latest'
HTMLlib = '/data/lv1/pcottle/HTML'
figurelib = '/data/lv1/WebData'


homepage = os.path.join(HTMLlib,'miniMPL-homepage.html')

try:
    sys.path.append(codelib)
    from MPLcode import MPLtools as mtools
    from MPLcode import MPL_plot as mplot   
except ImportError:
    
    try:
        import MPLtools as mtools
        import MPL_plot as mplot
    
    except ImportError:
        raise Exception('You havent specified where your modules are located!')

#get current date and time
current_time = datetime.datetime.now()

#get time tag of "latest.png" (set to current time if latest.png not found)
try:
    latestfig_time = os.path.getctime(figurelib+'latest.png')
except OSError:
    latestfig_time = current_time

#if it's been more than an hour since the last figure was generated, 
#generate a new one - otherwise do nothing    
if current_time-latestfig_time > datetime.timedelta(hours=1):
    mplfiles = []
    
    os.chdir(datalib)
    
    # gather all mpl files from the latest data folder and sort them by name
    #note: nnames are in the format "yyyymmddhh##.mpl" so this is also a de-facto
    #sort for creation date
    
    for temp in glob.glob("*.mpl"):
        mplfiles.append(temp)
    
    mplfiles.sort()
    
    # take files from last 24 hours, 
    
    mpl_selected = []
    
    for tempfile in reversed(mplfiles):
        temptime = os.path.getctime(tempfile)
        if current_time-temptime < datetime.timedelta(day=1):
            mpl_selected.append(tempfile)
        else:
            break
        

    
    #set altitude range and time step sizes
    
    altrange = np.arange(150,15000,30)#meters
    timestep = '60S' #seconds
    
    
    for f in mpl_selected:
        MPLdat_temp = mtools.MPL()
        MPLdat_temp.fromMPL(f)
        MPLdat_temp.alt_resample(altrange)
    
        try:
            MPLdat_event.append(MPLdat_temp)
        except NameError:
            MPLdat_event = MPLdat_temp
       
    #sort by index to make certain data is in order then set date ranges to match
    MPLdat_event.header.sort_index()
    
    for n in range(MPLdat_event.header['numchans'][0]):
        data = MPLdat_event.data[n]
        data = data.sort_index()
    
    MPLdat_event.time_resample(timestep)
    MPLdat_event.range_cor()    
    MPLdat_event.calculate_NRB(showplots = False)
    MPLdat_event.calculate_depolrat()
       
    #Now generate a new figure.
    
    copol = MPLdat_event.NRB[0]
    depol = MPLdat_event.depolrat[0]
    
    copol_min = 0.0
    copol_max = 1.0
    copol_step = 5
    copolticks = np.linspace(copol_min,copol_max,copol_step)
    
    depol_min = 0.0
    depol_max = 0.5
    depol_step = 5
    depolticks = np.linspace(depol_min,depol_max,depol_step)
    
    alt = copol.columns
    time_index = copol.index
    
    hours = ['03','06','09','12','15','18','21']
    
    fsize = 18 #baseline font size
    ar = 2.0  #aspect ratio
    figheight = 12 #inches
    
    
    plt.rc('font', family='serif', size=fsize)
    
    
    fig = plt.figure()
    
    h_set = range(1,25)
    h_set = map(str,h_set)
    
    ax1 = fig.add_subplot(2,1,1)
    im1 = mplot.backscatter_plot(fig, ax1, ar, time_index,alt[::-1],copol.T[::-1], (copol_min,copol_max), fsize = fsize)
    cbar1 = fig.colorbar(im1, orientation = 'vertical', aspect = 6)
    cbar1.set_ticks(copolticks)
    cbar1.set_ticklabels(copolticks)
    cbar1.ax.tick_params(labelsize = fsize)
    cbar1.ax.set_ylabel('$[counts*km^{2}/(\mu s*\mu J)$')
    mplot.dateticks(ax1, time_index, hours = hours, fsize = fsize, tcolor = 'w')
    ax1.set_xticklabels([])
    t1 = ax1.set_title('Normalized Relative Backscatter', fontsize = fsize+10)
    t1.set_y(1.03)
            
    ax2 = fig.add_subplot(2,1,2)
    im2 = mplot.depol_plot(fig, ax2, ar, time_index,alt[::-1],depol.T[::-1], (depol_min,depol_max), fsize = fsize)
    cbar2 = fig.colorbar(im2, orientation = 'vertical', aspect = 6)
    cbar2.set_ticks(depolticks)
    cbar2.set_ticklabels(depolticks)
    cbar2.ax.tick_params(labelsize = fsize)
    #set axis ranges and tickmarks based on data ranges
    mplot.dateticks(ax2, time_index, hours = hours, fsize = fsize)
    ax2.set_xlabel('Time [Local]',fontsize = fsize+4)
    fig.autofmt_xdate()
    t2 = ax2.set_title('Linear Depolarization Ratio',fontsize = fsize+10)
    t2.set_y(1.03)    
    
    
    fig.set_size_inches(figheight*ar,figheight) 
         
    #save fig as "latest.png"
    plt.savefig(figurelib+'latest.png',dpi = 100, edgecolor = 'b', bbox_inches = 'tight')
    
    
    #update html file
    with open(homepage) as html:
        for line in html:
            if "<h2>Latest Image:" in line:
                year = str(current_time.year)
                month = str(current_time.month)
                day = str(current_time.day)
                hour = str(current_time.hour)
                minute = str(current_time.minute)
                
                line = "<h2>Latest Image: "+year+"-"+month+"-"+day+" at "+hour+":"+minute+"</h2>"
                





