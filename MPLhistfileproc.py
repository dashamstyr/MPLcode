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
sns.set_context('poster')

if sys.platform == 'win32': 
    import matplotlib
    from matplotlib import pyplot as plt
    datadir='E:\CORALNet\ASCII_Files\Smoke2012\UBC\August\Processed'
    topdir='E:\CORALNet\ASCII_Files\Smoke2012\UBC\August\Processed'
    codelib = 'C:\Users\dashamstyr\Dropbox\Python_Scripts\GIT_Repos\MPLcode'
    figurelib = 'E:\CORALNet\ASCII_Files\Smoke2012\UBC\August\Figures'
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
timestep='600S'
SNRthresh=0.0
loadfiletype='LNCHDF'
savetype='HDF5'
saveresult=True
doplots=False
recalc=False


#print 'Performing initial file pre-precessing ...'
#mplfiles=mhist.filegetter(filedir=datadir,savedir=topdir,starttime=startdate,
#                          endtime=enddate,timestep=timestep,altitudes=altitudes,
#                          topickle=False,doprogress=True,recalc=False,verbose=False,
#                          SNRmask=False)
#mplfiles.sort()
#print ' ... Pre-processing Done!'
os.chdir(topdir)
mplfiles=mtools.get_files('Select LNC file',filetype = ('.h5','*.h5'))

filttype='Sub-Type'
#filtdict={'Dust':['Dust'],'Polluted_Dust':['Polluted Dust'],'Smoke_Urban':['Smoke / Urban'],
#          'Water_Cloud':['Water Cloud'],'Mixed_Cloud':['Mixed Cloud'],'Ice_Cloud':['Ice Cloud']}

filtdict={'Smoke':['Smoke / Urban','Polluted Dust']}
for layertype,filtlist in filtdict.iteritems():
    
#    datatypes={'NRB':['NRBhistdat_{0}'.format(layertype),(0.0,5.0),'NRBhistplot_{0}.png'.format(layertype)],
#                  'Depolrat':['Depolhistdat_{0}'.format(layertype),(0.0,1.0),'Depolhistplot_{0}.png'.format(layertype)],
#                  'SNR':['SNRhistdat_{0}'.format(layertype),(0.0,100.0),'SNRhistplot_{0}.png'.format(layertype)]}
#    datatypes={'Delta':['Deltahistdat_{0}'.format(layertype),(0.0,5.0),'Deltahistplot_{0}.png'.format(layertype)]}              
    datatypes={'BR':['BRhistdat_{0}'.format(layertype),(0.0,20.0)],'PR':['PRhistdat_{0}'.format(layertype)]}
    print 'Extracting {0} layer data'.format(layertype)
    for dtype,dspecs in datatypes.iteritems(): 
                
        print 'Generating {0} data for {1} layer type'.format(dtype,layertype)
        n=0
        fulldict=None
        dftotal=None
        
        for filename in mplfiles:
            mhist.progress(n,len(mplfiles))
            n+=1
            dffilt=mhist.dataextractor(datatype=filttype,loadfiletype=loadfiletype,loadfilename=filename,toHDF=False)
            dfSNR=mhist.dataextractor(datatype='SNR',loadfiletype=loadfiletype,loadfilename=filename,toHDF=False)
            df=mhist.dataextractor(datatype=dtype,loadfiletype=loadfiletype,loadfilename=filename,toHDF=False)
            dfplot=mhist.scenefilter(df,dffilt,filterterms=filtlist,filtmode='Type')
            dfplot=mhist.scenefilter(dfplot,dfSNR,filterterms=[SNRthresh],filtmode='GT')
            if savetype is 'HDF5':
                if dftotal is None:
                    dftotal=copy.deepcopy(dfplot)
                else:
                    dftotal=dftotal.append(dfplot)
            else:
                dftemp=dfplot.fillna(-99999,inplace=True)
                tempdict=mhist.histplot1D(dftemp,datatype='df',numbins=100,binrange=dspecs[1],
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
                store.close()
        if doplots:
            print 'Generating Histogram plots for {0} data on {1} layer type'.format(dtype,layertype)
            if fulldict is not None:                            
                fulldict=mhist.histplot1D(fulldict,datatype='histdict',numbins=100,doplot=True,saveplot=True,plotfilename=os.path.join(topdir,dspecs[2]))
            if dftotal is not None:
                histdat=dftotal.stack().dropna()
                fig=plt.figure()
                ax=fig.add_subplot(111)
                ax=sns.distplot(histdat,rug=False,kde=False,norm_hist=True)
                fig.canvas.print_figure(os.path.join(topdir,dspecs[2]),dpi = 100, edgecolor = 'b', bbox_inches = 'tight') 