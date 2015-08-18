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
import copy
import pickle 
import seaborn as sns

if sys.platform == 'win32': 
    import matplotlib
    from matplotlib import pyplot as plt
    datadir='C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Ucluelet Files\\testfile\\'
    topdir='C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Ucluelet Files\\testfile\Stats\\'
    codelib = 'C:\Users\dashamstyr\Dropbox\Python_Scripts\GIT_Repos\MPLcode'
    figurelib = 'K:\All_MPL_Data\\'
    systimezone = pytz.timezone('US/Pacific')
else:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    datadir='/data/lv1/pcottle/MPLData/Raw_Files'
    topdir='/data/lv1/pcottle/MPLStats/'
    codelib = '/data/lv1/pcottle/MPLCode/'
    systimezone = pytz.utc

localzone = pytz.timezone('US/Pacific')

try:
    sys.path.append(codelib)
    from MPLcode import MPLtools as mtools
    from MPLcode import MPLprocesstools as mproc
    from MPLcode import MPLhisttools as mhist
except ImportError:
    
    try:
        import MPLtools as mtools
        import MPLprocesstools as mproc
        import MPLhisttools as mhist
    except ImportError:
        raise Exception('You havent specified where your modules are located!')
   
startdate=datetime.datetime(2011,4,23,00)
enddate=datetime.datetime(2018,4,23,1)
altitudes=np.arange(0.150,10.000,0.030)
timestep='240S'
SNRthresh=2.0
dfcollect=True
saveresult=True
doplots=True

print 'Performing initial file pre-precessing ...'
mplfiles=mhist.filegetter(filedir=datadir,savedir=topdir,starttime=startdate,
                          endtime=enddate,timestep=timestep,altitudes=altitudes,
                          topickle=False,doprogress=True,recalc=False,verbose=False,
                          SNRmask=False)
mplfiles.sort()
print ' ... Pre-processing Done!'


filtlist=['Aerosol']

filtdict={'Aerosol':['Aerosol'],'Cloud':['Cloud'],'PBL':['PBL'],'All':['PBL','Cloud','Aerosol']}

for layertype,filtlist in filtdict.iteritems():
    
    datatypes={'NRB':['NRBhistdat_{0}'.format(layertype),(0.0,5.0),'NRBhistplot_{0}.png'.format(layertype)],
                  'Depolrat':['Depolhistdat_{0}'.format(layertype),(0,1.0),'Depolhistplot_{0}.png'.format(layertype)],
                  'SNR':['SNRhistdat_{0}'.format(layertype),(0,100.0),'SNRhistplot_{0}.png'.format(layertype)]}
                  
    print 'Extracting {0} layer data'.format(layertype)
    for dtype,dspecs in datatypes.iteritems(): 
        print 'Generating {0} data for {1} layer type'.format(dtype,layertype)
        n=0
        fulldict=None
        dftotal=None
        for filename in mplfiles:
            mhist.progress(n,len(mplfiles))
            n+=1
            dffilt=mhist.dataextractor(datatype='Type',loadfiletype='HDF5',loadfilename=filename,toHDF=False)
            dfSNR=mhist.dataextractor(datatype='SNR',loadfiletype='HDF5',loadfilename=filename,toHDF=False)
            df=mhist.dataextractor(datatype=dtype,loadfiletype='HDF5',loadfilename=filename,toHDF=False)
            dfplot=mhist.scenefilter(df,dffilt,filterterms=filtlist)
            dfplot=mhist.SNRfilter(df,dfSNR,thresh=SNRthresh)
            dfplot.fillna(value=-99999,inplace=True)
            if dfcollect:
                if dftotal is None:
                    dftotal=dfplot
                else:
                    dftotal.append(dfplot)
            else:
                tempdict=mhist.histplot1D(dfplot,datatype='df',numbins=100,binrange=dspecs[1],
                                            doplot=False,saveplot=False,cumulative=False,normalize=False)
                if fulldict is None:
                    fulldict=copy.deepcopy(tempdict)
                else:
                    fulldict['counts']+=tempdict['counts']
        
        print '{0} data generation for {1} layer type complete!'.format(dtype,layertype) 
        if saveresult:
            print 'Saving results to file {0}'.format(dspecs[0])
            if fulldict is not None:
                savefilename='{0}.p'.format(dspecs[0])
                pickle.dump(fulldict,open(os.path.join(topdir,savefilename),'wb'))
            if dftotal is not None:
                savefilename='{0}.h5'.format(dspecs[0])
                store=pan.HDFStore(os.path.join(topdir,savefilename))
                store['histdat']=dftotal
        if doplots:
            print 'Generating Histogram plots for {0} data on {1} layer type'.format(dtype,layertype)
            if fulldict is not None:                            
                fulldict=mhist.histplot1D(fulldict,datatype='histdict',numbins=100,doplot=True,saveplot=True,plotfilename=os.path.join(topdir,dspecs[2]))
            if dftotal is not None:
                fig=sns.distplot(df, numbins=100)
                fig.canvas.print_figure(os.path.join(topdir,dspecs[2]),dpi = 100, edgecolor = 'b', bbox_inches = 'tight') 