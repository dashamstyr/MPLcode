# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 17:04:59 2014
Updated on May 7,2014
@author: dashamstyr

Program to generate processed image of files from most recent 24 hour window in miniMPL data file.

"""

import os, sys
import pandas as pan
import datetime,time
import numpy as np
import glob


if sys.platform == 'win32': 
    import matplotlib
    from matplotlib import pyplot as plt
    datalib = 'C:\Users\dashamstyr\Dropbox\MPL website\pcottle\MPLData\Latest'
    HTMLlib = 'C:\Users\dashamstyr\Dropbox\MPL website\WebData'
    codelib = 'C:\Users\dashamstyr\Dropbox\MPL website\pcottle\MPLCode'
    figurelib = 'C:\Users\dashamstyr\Dropbox\MPL website\WebData'
    template = 'C:\Users\dashamstyr\Dropbox\MPL website\pcottle\HTML\index_temp.html'
else:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    codelib = '/data/lv1/pcottle/MPLCode/'
    datalib = '/data/lv1/pcottle/MPLData/Latest/'
    HTMLlib = '/data/lv1/WebData/'
    figurelib = '/data/lv1/WebData/'
    template = '/data/lv1/pcottle/HTML/index_temp.html'

homepage = os.path.join(HTMLlib,'index.html')
latest_figname = os.path.join(figurelib,'latest.png')

try:
    sys.path.append(codelib)
    from MPLcode import MPLtools as mtools
    from MPLcode import MPLplot as mplot   
except ImportError:
    
    try:
        import MPLtools as mtools
        import MPLplot as mplot
    
    except ImportError:
        raise Exception('You havent specified where your modules are located!')

#get current date and time
current_time = datetime.datetime.now()

#get time tag of "latest.png" (set to current time if latest.png not found)
if os.path.isfile(latest_figname):        
    latestfig_time = datetime.datetime.fromtimestamp(os.path.getctime(latest_figname))
else:
    latestfig_time = datetime.datetime.fromtimestamp(0)

#if it's been more than an hour since the last figure was generated, 
#generate a new one - otherwise do nothing    
if current_time-latestfig_time > datetime.timedelta(hours=1):
    mplfiles = []
    
    olddir = os.getcwd()
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
        temptime = datetime.datetime.fromtimestamp(os.path.getctime(tempfile))
        if current_time-temptime < datetime.timedelta(hours=24):
            mpl_selected.append(tempfile)
        else:
            continue
    
    if len(mpl_selected) > 0:
    #set altitude range and time step sizes
    
        altrange = np.arange(150,15000,30)#meters
        timestep = '60S' #seconds
        
        
        for f in mpl_selected:
            MPLdat_temp = mtools.MPL()
            MPLdat_temp.fromMPL(f)
            MPLdat_temp.alt_resample(altrange, verbose=False)
        
            try:
                MPLdat_event.append(MPLdat_temp)
            except NameError:
                MPLdat_event = MPLdat_temp
           
        #sort by index to make certain data is in order then set date ranges to match
        MPLdat_event.header.sort_index()
        
        for n in range(MPLdat_event.header['numchans'][0]):
            data = MPLdat_event.data[n]
            data = data.sort_index()
        
        MPLdat_event.time_resample(timestep,verbose=False)
        MPLdat_event.range_cor()    
        MPLdat_event.calculate_NRB(showplots = False)
        MPLdat_event.calculate_depolrat()
           
        #Now generate a new figure.       
        kwargs = {'saveplot':True,'showplot':False,'verbose':False,
                    'savefilepath':figurelib,'savefilename':latest_figname}
        
        mplot.doubleplot(MPLdat_event,**kwargs)

        
        #update html file
        with open(template,'r') as htmltemp:
            data = htmltemp.readlines()
        for n in range(len(data)):
            if "<h2>Latest Image:" in data[n]:
                year = str(current_time.year)
                month = str(current_time.month).zfill(2)
                day = str(current_time.day).zfill(2)
                hour = str(current_time.hour).zfill(2)
                minute = str(current_time.minute).zfill(2)
                
                data[n] = '<h2>Latest Image: {0}-{1}-{2} at {3}:{4}</h2>'.format(year,month,day,hour,minute)                
        with open(homepage,'w') as htmlhome:
            htmlhome.writelines(data)
            
    os.chdir(olddir)
    
    
    
    
