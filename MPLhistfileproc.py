# -*- coding: utf-8 -*-
"""
Created on Thu Oct 09 10:12:06 2014

@author: dashamstyr
"""

import os, sys
import pandas as pan
import datetime,time,pytz
import numpy as np
import glob
import MPLhisttools as mhist 

if sys.platform == 'win32': 
    import matplotlib
    from matplotlib import pyplot as plt
    topdir='C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA'
    savefile='C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\histpickle_all.p'
    HDFfile='C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\histdat_all.h5'
    datalib = 'C:\Users\dashamstyr\Dropbox\MPL website\pcottle\MPLData\Latest'
    sitelib = 'C:\Users\dashamstyr\Dropbox\MPL website\WebData'
    codelib = 'C:\Users\dashamstyr\Dropbox\Python_Scripts\GIT_Repos\MPLcode'
    figurelib = 'C:\Users\dashamstyr\Dropbox\MPL website\WebData'
    systimezone = pytz.timezone('US/Pacific')
else:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    topdir='/data/lv1/sigma/'
    savefile='/data/lv1/pcottle/MPLData/histpickle_Sigma.p'
    HDFfile='/data/lv1/pcottle/MPLData/histdat_Sigma.h5'
    codelib = '/data/lv1/pcottle/MPLCode/'
    datalib = '/data/lv1/pcottle/MPLData/Latest/'
    sitelib = '/data/lv1/WebData/'
    figurelib = '/data/lv1/WebData/'
    systimezone = pytz.utc

localzone = pytz.timezone('US/Pacific')

#try:
#    sys.path.append(codelib)
#    from MPLcode import MPLtools as mtools
#    from MPLcode import MPLprocesstools as mproc
#    from MPLcode import MPLhisttools as mhist
#except ImportError:
    
#    try:
#        import MPLtools as mtools
#        import MPLprocesstools as mproc
#        import MPLhisttools as mhist
#    except ImportError:
#        raise Exception('You havent specified where your modules are located!')
   
startdate=datetime.datetime(2013,4,23,00)
enddate=datetime.datetime(2015,4,23,1)
altitudes=np.arange(150,15000,30)
timestep='600S'
mplfiles=mhist.filegetter(filedir=topdir,altitudes=altitudes,starttime=startdate,
                          endtime=enddate,timestep=timestep,topickle=True,picklefile=savefile,
                          verbose=True)

#with open(savefile,'rb') as sfile:
#    allfiles=pickle.load(sfile)
#
#datatypes={'NRB':['NRBhistdat.h5',(0.0,5.0)],'Depolrat':['Depolhistdat.h5',(0,1.0)],'SNR':['SNRhistdat.h5',(0,10.0)]}
#alldict={}
#n=0
#for dtype,dspecs in datatypes.iteritems():
#    df=dataextractor(datatype=dtype,fromHDF=True,loadfilename=os.path.join(topdir,dspecs[0]),toHDF=False)
#    dfplot=df[df>0]
#    dfplot.fillna(value=-99999,inplace=True)    
#    
#    alldict[dtype]=mhist.histplot1D(df,numbins=100,binrange=dspecs[1],xlabel=dtype,fignum=n,cumulative=False,normalize=True)
#    n+=1