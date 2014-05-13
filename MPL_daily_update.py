# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 13:04:21 2014
Updated on May 7,2014
@author: User
"""
import os, sys
import pandas as pan
import datetime,time
import numpy as np
import glob



#Step 1: Set folders for locations of unsorted raw data files, python code, webdata folder
if sys.platform == 'win32':   
    import matplotlib
    from matplotlib import pyplot as plt    
    datalib = 'C:\Users\dashamstyr\Dropbox\MPL website\pcottle\MPLData\Raw_Files'
    sitelib = 'C:\Users\dashamstyr\Dropbox\MPL website\WebData'
    codelib = 'C:\Users\dashamstyr\Dropbox\MPL website\pcottle\MPLCode'
    historyfile = 'C:\Users\dashamstyr\Dropbox\MPL website\pcottle\MPLData\MPL_location_history.txt'
else:    
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    datalib = '/data/lv1/pcottle/MPLData/Raw_Files'
    sitelib = '/data/lv1/WebData'
    codelib = '/data/lv1/pcottle/MPLCode'
    historyfile = '/data/lv1/pcottle/MPLData/MPL_location_history.txt'

try:
    sys.path.append(codelib)
    from MPLcode import MPLtools as mtools
    reload(MPLtools)
    from MPLcode import HTMLtools as htools 
    reload(HTML_tools)
    from MPLcode import MPLplot as mplot
except ImportError:
    
    try:
        import MPLtools as mtools
        reload(mtools)
        import HTMLtools as htools
        reload(htools)
        import MPLplot as mplot
    
    except ImportError:
        raise Exception('You havent specified where your modules are located!')

olddir=os.getcwd()
#Set global variables for plots and data files
#--------------------------------------------------------------------
#set altitude range and time step sizes for processing and plots

altrange = np.arange(150,15000,30)#meters
timestep = '60S' #seconds

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

for f,date in zip(filenames,filedates):
    
    #Step 4a: Bundle files by day and generate new processed filename YYYYMMDD.h5 (&.png) 
    
    
    if date == newday:
        dailympl.append(os.path.join(datalib,f))
    else:
        tempyear = date[0]
        tempmonth = date[1]    
        tempday = date[2]
        #Check if .h5 file exists
        
        savefilename = '{0}{1}{2}'.format(tempyear,tempmonth,tempday)
        
        if savefilename in procfiles:
            continue
        else:        
            saveloc = htools.findloc(f,location_dict)
            savefoldername = '{0}-{1}'.format(tempyear,tempmonth)
            
            loclib = os.path.join(sitelib,saveloc)
            procdatlib = os.path.join(loclib,'Data_Files',savefoldername)
            procimglib = os.path.join(loclib,'Image_Files',savefoldername)           
            
            if not os.path.isdir(procdatlib):
                os.mkdir(procdatlib)
            if not os.path.isdir(procimglib):
                os.mkdir(procimglib)
#            MPLdat_event = mtools.MPL()
            
            for mpl in dailympl:
                MPLdat_temp = mtools.MPL()
                MPLdat_temp.fromMPL(mpl)
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
            'savefilepath':procimglib,'savefilename':'{0}.png'.format(savefilename)}
        
            mplot.doubleplot(MPLdat_event,**kwargs)
            
            plt.close('all')

            #Step 4c: Switch to YYYY-MM folder (create if not there already) and save .h5, then repeat for .png
            
            os.chdir(procdatlib)
            MPLdat_event.save_to_HDF('{0}.h5'.format(savefilename))
            
            newday = date    
            dailympl = []
            dailympl.append(os.path.join(datalib,f))
            os.chdir(olddir)


