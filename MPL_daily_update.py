# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 13:04:21 2014

@author: User
"""
import os, sys
import pandas as pan
import datetime,time
import numpy as np
import glob
import matplotlib
from matplotlib import pyplot as plt


#Step 1: Set folders for locations of unsorted raw data files, python code, webdata folder
if sys.platform == 'win32':    
#    matplotlib.use('Agg')    
    datalib = 'C:\Users\dashamstyr\Dropbox\MPL website\pcottle\MPLData\Raw_Files'
    sitelib = 'C:\Users\dashamstyr\Dropbox\MPL website\WebData'
    codelib = 'C:\Users\dashamstyr\Dropbox\MPL website\pcottle\MPLCode'
    historyfile = 'C:\Users\dashamstyr\Dropbox\MPL website\pcottle\MPLData\MPL_location_history.txt'
else:
    matplotlib.use('Agg')
    datalib = '/data/lv1/pcottle/MPLData/Raw_Files'
    sitelib = '/data/lv1/WebData'
    codelib = '/data/lv1/pcottle/MPLCode'
    historyfile = '/data/lv1/pcottle/MPLData/MPL_location_history.txt'

try:
    sys.path.append(codelib)
    from MPLcode import MPLtools as mtools
    reload(MPLtools)
    from MPLcode import HTML_tools as htools 
    reload(HTML_tools)
    from MPLcode import MPL_plot as mplot
except ImportError:
    
    try:
        import MPLtools as mtools
        reload(mtools)
        import HTML_tools as htools
        reload(htools)
        import MPL_plot as mplot
    
    except ImportError:
        raise Exception('You havent specified where your modules are located!')

olddir=os.getcwd()
#Set global variables for plots and data files
#--------------------------------------------------------------------
#set altitude range and time step sizes for processing and plots

altrange = np.arange(150,15000,30)#meters
timestep = '60S' #seconds

#set ranges and tickmarks for plots
copol_min = 0.0
copol_max = 1.0
copol_step = 5
copolticks = np.linspace(copol_min,copol_max,copol_step)

depol_min = 0.0
depol_max = 0.5
depol_step = 5
depolticks = np.linspace(depol_min,depol_max,depol_step)

# set hour markers for time axis
hours = ['03','06','09','12','15','18','21']
h_set = range(1,25)
h_set = map(str,h_set)


#set font sizes, styles and figure size and aspect ratio
        
fsize = 18 #baseline font size
ar = 2.0  #aspect ratio
figheight = 12 #inches

plt.rc('font', family='serif', size=fsize)

#----------------------------------------------------------------------


location_dict = htools.locations(historyfile)

for loc,dates in location_dict.iteritems():
    templib = os.path.join(sitelib,loc)
    tempdatlib = os.path.join(templib,'Data_Files')
    tempimglib = os.path.join(templib,'Image_Files')
    
    for lib in [templib,tempdatlib,tempimglib]:
        if not os.path.isdir(lib):
            os.mkdir(lib)
        
    if 'na' in dates[1]:
        current_loc = loc
        
        
   
#Step 2: Open raw data folder and make a list of all .mpl files
#os.chdir(datalib)
allfiles = os.listdir(datalib)


procfiles=[]
for root,dirs,files in os.walk(sitelib):
    for f in files:
        if f.endswith('.h5'):
            procfiles.append(f.split('.')[0])
            
#Step 3: Make dict of all year-month-day combos in .mpl files
rawfiles=[]
for f in allfiles:
    if f.endswith('.mpl') and f[:8] not in procfiles:
        rawfiles.append(f)


filenames,filedates = htools.listdates(rawfiles)
    
#Step 4: Loop through .mpl files and ...

newday = filedates[0]
dailympl =[]

figcount=0
for f,date in zip(filenames,filedates):
    
    #Step 4a: Bundle files by day and generate new processed filename YYYYMMDD.h5 (&.png) 
    
    
    if date == newday:
        dailympl.append(os.path.join(datalib,f))
    else:
        tempyear = date[0]
        tempmonth = date[1]    
        tempday = date[2]
        #Check if .h5 file exists
        
        savefilename = tempyear+tempmonth+tempday
        
        if savefilename in procfiles:
            continue
        else:        
            saveloc = htools.findloc(f,location_dict)
            savefoldername = tempyear+'-'+tempmonth
            
            loclib = os.path.join(sitelib,saveloc)
            procdatlib = os.path.join(loclib,'Data_Files',savefoldername)
            procimglib = os.path.join(loclib,'Image_Files',savefoldername)           
            
            if not os.path.isdir(procdatlib):
                os.mkdir(procdatlib)
            if not os.path.isdir(procimglib):
                os.mkdir(procimglib)
#            MPLdat_event = mtools.MPL()
            
            for mpl in dailympl:
                print os.getcwd()
                MPLdat_temp = mtools.MPL()
                MPLdat_temp.fromMPL(mpl)
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
            
            alt = copol.columns
            time_index = copol.index
            
            fig = plt.figure(figcount)        
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
            #Step 4c: Switch to YYYY-MM folder (create if not there already) and save .h5, then repeat for .png
            
            os.chdir(procdatlib)
            MPLdat_event.save_to_HDF(savefilename+'.h5')
            os.chdir(procimglib)
            fig.canvas.print_figure(savefilename+'.png',dpi = 100, edgecolor = 'b', bbox_inches = 'tight')
            
            #reset newday, dailympl,and MPLdat_event for next iteration
            fig.clf()
            del fig, MPLdat_event
            newday = date    
            dailympl = []
            dailympl.append(os.path.join(datalib,f))
            figcount+=1
            os.chdir(olddir)


