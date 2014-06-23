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
    codelib = 'C:\Users\dashamstyr\Dropbox\Python_Scripts\GIT_Repos\MPLcode'
    historyfile = 'C:\Users\dashamstyr\Dropbox\MPL website\pcottle\MPLData\MPL_location_history.txt'
else:    
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    datalib = '/data/lv1/pcottle/MPLData/Raw_Files'
    sitelib = '/data/lv1/WebData'
    codelib = '/data/lv1/pcottle/MPLCode'
    historyfile = '/data/lv1/pcottle/MPLData/MPL_location_history.txt'

homepage = os.path.join(sitelib,'index.html')

try:
    sys.path.append(codelib)
    from MPLcode import MPLtools as mtools
    from MPLcode import HTMLtools as htools 
    from MPLcode import MPLplot as mplot
    from MPLcode import MPLprocesstools as mproc
except ImportError:
    
    try:
        import MPLtools as mtools
        import HTMLtools as htools
        import MPLplot as mplot
        import MPLprocesstools as mproc    
    except ImportError:
        raise Exception('You havent specified where your modules are located!')

olddir=os.getcwd()
#Set global variables for plots and data files
#--------------------------------------------------------------------
#set altitude range and time step sizes for processing and plots

altrange = np.arange(150,15060,60)#meters
timestep = '120S' #seconds

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

if rawfiles:
    np.sort(rawfiles)
    filenames,filedates = htools.listdates(rawfiles)
        
    #Step 4: Loop through .mpl files and ...
    
    newday = filedates[0]
    dailympl =[]
    MPLdat_event = []
    
    for f,date in zip(filenames,filedates):
        
        #Step 4a: Bundle files by day and generate new processed filename YYYYMMDD.h5 (&.png) 
        if date != newday or f==filenames[-1]:
            tempyear = newday[0]
            tempmonth = newday[1]    
            tempday = newday[2]
            
            savefilename = '{0}{1}{2}'.format(tempyear,tempmonth,tempday)
            
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
            
                if MPLdat_event:
                    MPLdat_event.append(MPLdat_temp)
                else:
                    MPLdat_event = MPLdat_temp
               
            #sort by index to make certain data is in order then set date ranges to match
            MPLdat_event.header.sort_index()
            
            for n in range(MPLdat_event.header['numchans'][0]):
                data = MPLdat_event.data[n]
                data = data.sort_index()
            
            MPLdat_event.time_resample(timestep,verbose=False)
            MPLdat_event.calc_all()
            
            MPLdat_event=mproc.NRB_mask_all(MPLdat_event)
                
            
               
            #Now generate a new figure.
               
            kwargs = {'saveplot':True,'showplot':False,'verbose':False,
            'savefilepath':procimglib,'savefilename':'{0}.png'.format(savefilename),
            'SNRmask':True}

            mplot.doubleplot(MPLdat_event,**kwargs)
            
            plt.close('all')
    
            #Step 4c: Switch to YYYY-MM folder (create if not there already) and save .h5, then repeat for .png
            
            os.chdir(procdatlib)
            MPLdat_event.save_to_HDF('{0}.h5'.format(savefilename))
            
            newday = date    
            dailympl = []
            MPLdat_event = []
            dailympl.append(os.path.join(datalib,f))
            os.chdir(olddir)
        else:
            dailympl.append(os.path.join(datalib,f))
            newday = date
    
    
    with open(homepage,'r') as htmltemp:
        data = htmltemp.readlines()
    for n in range(len(data)):
        if "<h1>Current Location:" in data[n]:
            data[n] = '<h1>Current Location: {0}</h1>\r\n'.format(current_loc)
    with open(homepage,'w') as htmlhome:
        htmlhome.writelines(data)
