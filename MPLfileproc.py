#from __future__ import absolute_import
import os,glob,site,sys

import numpy as np
import datetime
import pandas as pan
from scipy import signal

import MPLtools as mtools
import MPLprocesstools as mproc
import MPLplot as mplot

if sys.platform == 'win32': 
    from matplotlib import pyplot as plt
    datalib = 'C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Whistler-0327'
    proclib = 'C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Whistler-0327\Processed'
    figurelib = 'C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Whistler-0327\Figures'
else:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    lib = '/data/lv1/pcottle'
    datalib = '/data/lv1/pcottle/MPLData/Smoke2015'
    proclib = '/data/lv1/pcottle/MPLData/Smoke2015/Processed'
    figurelib = '/data/lv1/pcottle/MPLData/Smoke2015/Figures'

def fileproc(**kwargs):

    """
    ----------------------------------------------------------------------------
    Uses tools created in MPL_tools to open all files in a folder and resample
    them to a regular spacing in altitude/date the concatenates them into one
    pandas dataframe and saves it as a 
    created: July 05, 2012
    updated: May 09,2014
    ----------------------------------------------------------------------------     
    inputs:
    **kwargs:
         altrange = a list of altitudes to process(becomes pandas columns)
         starttime = time point for first profile
         endtime = time point for last profile
         timestep = time between profiles
         interactive = boolean determining if user will select files or if all in folder will be processed
         doplot = boolean determining if plot will be made from processed file
         dolayers = boolean determining if the 
         saveplot = boolean determining if plot will be saved
         verbose = boolean determining if messages will be displayed
         
    
    """
    #general kwargs
    recalc=kwargs.get('recalc',True)
    rawfiles = kwargs.get('rawfiles',None)
    altrange = kwargs.get('altrange',np.arange(150,15000,30))        
    timestep = kwargs.get('timestep','60S')
    saveproc = kwargs.get('saveproc',True)
    doplot = kwargs.get('doplot',False)
    dolayers = kwargs.get('dolayers',False)
    docorrection = kwargs.get('docorrection',False)
    dolayerplot = kwargs.get('dolayerplot',False)
    docorplot = kwargs.get('docorplot',False)
    saveplot = kwargs.get('saveplot',False)
    savetype = kwargs.get('savetype','standard')
    showplot = kwargs.get('showplot',False)
    procsavepath = kwargs.get('procsavepath','.\Processed')
    plotsavepath = kwargs.get('plotsavepath','.\Figures')
    verbose = kwargs.get('verbose',False)
    #NRBmask kwargs
    NRBmask = kwargs.get('NRBmask',True)
    NRBmasktype = kwargs.get('NRBmasktype','profile')
    NRBthresh=kwargs.get('NRBthresh',3.0)
    NRBmin=kwargs.get('NRBmin',0.5)
    NRBminalt=kwargs.get('NRBminalt',0.150)
    NRBnumprofs=kwargs.get('NRBnumprofs',1)
    NRBwinsize=kwargs.get('NRBwinsize',3)
    NRBinplace=kwargs.get('NRBinplace',False)
    
    #findalllayers kwargs
    molthresh=kwargs.get('molthresh',1.0)
    noisethresh=kwargs.get('noisethresh',1.0) 
    depolnoisethresh=kwargs.get('depolnoisethresh',1.0)     
    bg_alt=kwargs.get('bg_alt',None)
    datatype=kwargs.get('datatype','NRB')
    winsize=kwargs.get('winsize',5)
    wavelet=kwargs.get('wavelet',signal.ricker)
    cloudthresh=kwargs.get('cloudthresh',(1.0,0.20))
    CWTwidth=kwargs.get('CWTwidth',2)
    minwidth=kwargs.get('minwidth',4)
    layerCWTrange=kwargs.get('layerCWTrange',np.arange(2,5))
    doPBL=kwargs.get('doPBL',True)
    PBLwavelet=kwargs.get('PBLwavelet',mproc.dog)
    PBLCWTrange=kwargs.get('PBLCWTrange',np.arange(2,10))
    PBLwidth=kwargs.get('PBLwidth',5)
    PBLnumprofs=kwargs.get('PBLnumprofs',1)
    PBLwinsize=kwargs.get('PBLwinsize',5)
    sigma0=kwargs.get('sigma0',None)
    depolsigma0=kwargs.get('depolsigma0',None)
    waterthresh=kwargs.get('waterthresh',0.10)
    icethresh=kwargs.get('icethresh',0.25)
    smokethresh=kwargs.get('smokethresh',0.05)
    dustthresh=kwargs.get('dustthresh',0.15)
    maxaeroalt=kwargs.get('maxaeroalt',10.0)
    
    #scene kwargs 
    PBLrat=kwargs.get('PBLrat',30.0)
    molrat=kwargs.get('molrat',0.0)
    moldepol=kwargs.get('moldepol',0.0035)
    
    #correction kwargs
    refalt=kwargs.get('refalt',None)
    calrange=kwargs.get('calrange',None)
    method=kwargs.get('method','klett2')
    lrat=kwargs.get('lrat',None)
    mode=kwargs.get('mode','copol')
    cor_numprofs=kwargs.get('cor_numprofs',1)
    cor_winsize=kwargs.get('cor_winsize',1)
    
    #plot kwargs
    NRB_limits = kwargs.get('NRB_limits',(0.0,1.0,0.2))  
    depol_limits = kwargs.get('depol_limits',(0.0,0.5,0.1))
    back_limits = kwargs.get('back_limits',(0.0,1e-7,2e-8))
    ext_limits = kwargs.get('ext_limits',(0.0,2e-6,4e-7))
    hours = kwargs.get('hours',['06','12','18']) 
    fsize = kwargs.get('fsize',18) #baseline font size
    SNRmask = kwargs.get('SNRmask',False)
    SNRthresh=kwargs.get('SNRthresh',1.0)
    SNRtype=kwargs.get('SNRtype','NRB')
    interpolate=kwargs.get('interpolate',None) #other methods include none, bilinear, gaussian and hermite
    #starttime and endtime are defined later      
    olddir=os.getcwd()
    
    if recalc:
        if rawfiles is None:
            rawfiles = mtools.get_files('Select files to process',filetype = ('.mpl','*.mpl'))
        
        [path,startfile] = os.path.split(rawfiles[0])
        [path,endfile] = os.path.split(rawfiles[-1])
        starttime = kwargs.get('starttime',MPLtodatetime(startfile))
        endtime = kwargs.get('endtime',MPLtodatetime(endfile))
            
        for r in rawfiles:
            [path,tempname] = os.path.split(r)
            if starttime <= MPLtodatetime(tempname) <= endtime: 
                MPLdat_temp = mtools.MPL()
                MPLdat_temp.fromMPL(r)
                MPLdat_temp.alt_resample(altrange,verbose=verbose)    
                try:
                    MPLdat_event.append(MPLdat_temp)
                except NameError:
                    MPLdat_event = MPLdat_temp
           
        #sort by index to make certain data is in order then set date ranges to match
        MPLdat_event.header.sort_index()
        
        for n in range(MPLdat_event.header['numchans'][0]):
            data = MPLdat_event.data[n]
            data = data.sort_index()
        
        MPLdat_event.time_resample(timestep,verbose=verbose)
        MPLdat_event.calc_all()
       
        if dolayers:
            layerkwargs={'timestep':timestep,'bg_alt':bg_alt,'datatype':datatype,
                         'molthresh':molthresh,'winsize':winsize,'noisethresh':noisethresh,
                         'depolnoisethresh':depolnoisethresh,'wavelet':wavelet,'cloudthresh':cloudthresh,
                         'CWTwidth':CWTwidth,'minwidth':minwidth,'layerCWTrange':layerCWTrange,'doPBL':doPBL,
                         'PBLwavelet':PBLwavelet,'PBLCWTrange':PBLCWTrange,'PBLwidth':PBLwidth,
                         'PBLnumprofs':PBLnumprofs,'PBLwinsize':PBLwinsize,'sigma0':sigma0,
                         'waterthresh':waterthresh,'icethresh':icethresh,'smokethresh':smokethresh,
                         'dustthresh':dustthresh,'sigma0':sigma0,'depolsigma0':depolsigma0,'maxaeroalt':maxaeroalt}
            layerdict=mproc.findalllayers(mplin=MPLdat_event,**layerkwargs)
            
            scenekwargs={'PBLrat':PBLrat,'molrat':molrat,'moldepol':moldepol,'cloudthresh':cloudthresh,
                         'waterthresh':waterthresh,'icethresh':icethresh,'smokethresh':smokethresh,
                         'dustthresh':dustthresh,'maxaeroalt':maxaeroalt}
            MPLdat_event=mproc.scenemaker(layerdict,**scenekwargs)
        
        if docorrection:
            corkwargs={'refalt':refalt,'calrange':calrange,'method':method,'lrat':lrat,'mode':mode,
                       'numprofs':cor_numprofs,'winsize':cor_winsize}
            
            MPLdat_event=mproc.basiccorrection(MPLdat_event,**corkwargs)
        
        if saveproc:
            if savetype=='standard':
                if len(rawfiles) == 1:   
                    d_filename = '{0}_proc.h5'.format(startfile.split('.')[0])
                else:        
                    d_filename = '{0}-{1}_proc.h5'.format(startfile.split('.')[0],endfile.split('.')[0])
                savepath=os.path.join(procsavepath,d_filename)
                try:
                    MPLdat_event.save_to_HDF(savepath)
                except IOError:
                    os.mkdir(procsavepath)
                    MPLdat_event.save_to_HDF(savepath)
            elif savetype=='IDL':
                if len(rawfiles) == 1:   
                    d_filename = '{0}_IDL.h5'.format(startfile.split('.')[0])
                else:        
                    d_filename = '{0}-{1}_IDL.h5'.format(startfile.split('.')[0],endfile.split('.')[0])
                savepath=os.path.join(procsavepath,d_filename)
                try:
                    MPLdat_event.save_to_IDL(procsavepath)
                except IOError:
                    os.mkdir(savepath)
                    MPLdat_event.save_to_IDL(savepath)
    else:
        if rawfiles is None:
            rawfiles = mtools.get_files('Select files to process',filetype = ('.h5','*.h5'))
        d_filename=rawfiles[0]
        MPLdat_event=mtools.MPL()
        MPLdat_event.fromHDF(d_filename)
        MPLdat_event.alt_resample(altrange,verbose=verbose) 
        
        
    if NRBmask:
        NRBmaskkwargs={'NRBmasktype':NRBmasktype,'NRBthreshold':NRBthresh,'NRBmin':NRBmin,
                       'minalt':NRBminalt,'numprofs':NRBnumprofs,'winsize':NRBwinsize,
                       'inplace':NRBinplace}
        MPLdat_event=mtools.NRB_mask_all(MPLdat_event,**NRBmaskkwargs)


    if doplot:
        plotfilename = '{0}.png'.format(d_filename.split('.')[0])
        plotkwargs={'altrange':altrange,'topplot_limits':NRB_limits,'bottomplot_limits':depol_limits,
                    'hours':hours,'fsize':fsize,'savefilename':plotfilename,'savefilepath':plotsavepath,
                    'SNRmask':SNRmask,'SNRthresh':SNRthresh,'SNRtype':SNRtype,'interpolate':interpolate,
                    'showplot':showplot,'saveplot':saveplot}
        
        mplot.doubleplot(MPLdat_event,plotfilename=plotfilename,**plotkwargs)
    
    if docorplot:
        coefplotfilename = '{0}-coefficients.png'.format(d_filename.split('.')[0])
        coefplotkwargs={'altrange':altrange,'toptype':'backscatter','bottomtype':'extinction',
                        'topplot_limits':back_limits,'bottomplot_limits':ext_limits,
                        'hours':hours,'fsize':fsize,'savefilename':coefplotfilename,
                        'savefilepath':plotsavepath,#'SNRmask':SNRmask,'SNRthresh':SNRthresh,
                        'interpolate':interpolate,'SNRtype':SNRtype,'showplot':showplot,'saveplot':saveplot}
            
        mplot.doubleplot(MPLdat_event,plotfilename=coefplotfilename,**coefplotkwargs)

    if dolayerplot:
        layerplotfilename = '{0}-layers.png'.format(d_filename.split('.')[0])
        layerplotkwargs={'altrange':altrange,'hours':hours,'fontsize':fsize,'plotfilename':layerplotfilename,
                         #'SNRmask':SNRmask,'SNRthresh':SNRthresh,'SNRtype':SNRtype,
                         'plotfilepath':plotsavepath,'showplot':showplot,'saveplot':saveplot}        
        mplot.colormask_plot(MPLdat_event,**layerplotkwargs)
#    if os.path.isdir(procsavepath):
#        os.chdir(procsavepath)
#    else:
#        os.makedirs(procsavepath)
#        os.chdir(procsavepath)
#    
#    if verbose:
#        print 'Saving '+d_filename
#    
#    if savetype=='standard':
#        MPLdat_event.save_to_HDF(d_filename)
#    elif savetype=='IDL':
#        MPLdat_event.save_to_IDL(d_filename)
    
    if verbose:
        print 'Done'
    
    os.chdir(olddir)
    return MPLdat_event
    
def getmplfiles(newdir,interactive=True, verbose=False):
    olddir = os.getcwd()
    os.chdir(newdir)
    
    if interactive:
        rawfiles = mtools.get_files('Select MPL files', filetype = ('.mpl', '*.mpl'))
    else:
        rawfiles = glob.glob('*.mpl')
    
    rawfiles.sort()
    
    if verbose:
        print rawfiles
    os.chdir(olddir)
    
    return rawfiles

def processall(**kwargs):
    newdir=kwargs.get('newdir',None)
    chunksize=kwargs.get('chunksize','days')
    startdate=kwargs.get('startdate',datetime.datetime(2009,01,01))
    enddate=kwargs.get('enddate',datetime.datetime.now())
    interactive=kwargs.get('interactive',False)
    #general kwargs
    altrange = kwargs.get('altrange',np.arange(150,15000,30))        
    timestep = kwargs.get('timestep','60S')
    doplot = kwargs.get('doplot',False)
    dolayers = kwargs.get('dolayers',False)
    docorrection = kwargs.get('docorrection',False)
    dolayerplot = kwargs.get('dolayerplot',False)
    docorplot = kwargs.get('docorplot',False)
    saveplot = kwargs.get('saveplot',False)
    savetype = kwargs.get('savetype','standard')
    showplot = kwargs.get('showplot',False)
    procsavepath = kwargs.get('procsavepath','.\Processed')
    plotsavepath = kwargs.get('plotsavepath','.\Figures')
    verbose = kwargs.get('verbose',False)
    #NRBmask kwargs
    NRBmask = kwargs.get('NRBmask',True)
    NRBmasktype = kwargs.get('NRBmasktype','profile')
    NRBthresh=kwargs.get('NRBthresh',3.0)
    NRBmin=kwargs.get('NRBmin',0.5)
    NRBminalt=kwargs.get('NRBminalt',0.150)
    NRBnumprofs=kwargs.get('NRBnumprofs',1)
    NRBwinsize=kwargs.get('NRBwinsize',3)
    NRBinplace=kwargs.get('NRBinplace',False)

    #findalllayers kwargs
    molthresh=kwargs.get('molthresh',1.0)
    noisethresh=kwargs.get('noisethresh',1.0)      
    depolnoisethresh=kwargs.get('depolnoisethresh',1.0)
    bg_alt=kwargs.get('bg_alt',None)
    datatype=kwargs.get('datatype','NRB')
    winsize=kwargs.get('winsize',5)
    wavelet=kwargs.get('wavelet',signal.ricker)
    noisethresh=kwargs.get('noisethresh',0.4)
    cloudthresh=kwargs.get('cloudthresh',(1.0,0.20))
    CWTwidth=kwargs.get('CWTwidth',2)
    minwidth=kwargs.get('minwidth',4)
    layerCWTrange=kwargs.get('layerCWTrange',np.arange(2,5))
    doPBL=kwargs.get('doPBL',True)
    PBLwavelet=kwargs.get('PBLwavelet',mproc.dog)
    PBLCWTrange=kwargs.get('PBLCWTrange',np.arange(2,10))
    PBLwidth=kwargs.get('PBLwidth',5)
    PBLnumprofs=kwargs.get('PBLnumprofs',1)
    PBLwinsize=kwargs.get('PBLwinsize',5)
    sigma0=kwargs.get('sigma0',None)
    depolsigma0=kwargs.get('depolsigma0',None)
    waterthresh=kwargs.get('waterthresh',0.10)
    icethresh=kwargs.get('icethresh',0.25)
    smokethresh=kwargs.get('smokethresh',0.05)
    dustthresh=kwargs.get('dustthresh',0.15)
    maxaeroalt=kwargs.get('maxaeroalt',10.0)
    
    #scene kwargs
    PBLrat=kwargs.get('PBLrat',30.0)
    molrat=kwargs.get('molrat',0.0)
    moldepol=kwargs.get('moldepol',0.0035)    
    
    #correction kwargs
    refalt=kwargs.get('refalt',None)
    calrange=kwargs.get('calrange',None)
    method=kwargs.get('method','klett2')
    lrat=kwargs.get('lrat',None)
    mode=kwargs.get('mode','copol')
    cor_numprofs=kwargs.get('cor_numprofs',1)
    cor_winsize=kwargs.get('cor_winsize',1)
    
    #plot kwargs
    NRB_limits = kwargs.get('NRB_limits',(0.0,1.0,0.2))  
    depol_limits = kwargs.get('depol_limits',(0.0,0.5,0.1))
    back_limits = kwargs.get('back_limits',(0.0,1e-7,2e-8))
    ext_limits = kwargs.get('ext_limits',(0.0,2e-6,4e-7))
    hours = kwargs.get('hours',['06','12','18']) 
    fsize = kwargs.get('fsize',18) #baseline font size
    SNRmask = kwargs.get('SNRmask',False)
    SNRthresh=kwargs.get('SNRthresh',1.0)
    SNRtype=kwargs.get('SNRtype','NRB')
    interpolate=kwargs.get('interpolate',None) #other methods include none, bilinear, gaussian and hermite
    
    olddir = os.getcwd()
    if not showplot:
        plt.ioff()
        
    if newdir is None:
        newdir = mtools.set_dir('Select directory to process')
    os.chdir(newdir)
    if interactive:
        rawfiles = mtools.get_files('Select MPL files', filetype = ('.mpl', '*.mpl'))
    else:
        rawfiles = glob.glob('*.mpl')
    
    rawfiles.sort()
    olddate=MPLtodatetime(rawfiles[0])
    tempfilelist=[rawfiles[0]]
    MPLlist=[]
    
    prockwargs={'altrange':altrange,'timestep':timestep,'doplot':doplot,'dolayers':dolayers,
                'docorrection':docorrection,'cor_numprofs':cor_numprofs,'cor_winsize':cor_winsize,
                'dolayerplot':dolayerplot,'docorplot':docorplot,'saveplot':saveplot,
                'savetype':savetype,'showplot':showplot,'procsavepath':procsavepath,
                'plotsavepath':plotsavepath,'verbose':verbose,'NRBmask':NRBmask,'NRBmasktype':NRBmasktype,
                'NRBthresh':NRBthresh,'NRBmin':NRBmin,'NRBminalt':NRBminalt,'NRBnumprofs':NRBnumprofs,'NRBwinsize':NRBwinsize,
                'NRBinplace':NRBinplace,'SNRmask':SNRmask,'SNRthresh':SNRthresh,'SNRtype':SNRtype,'molthresh':molthresh,
                'bg_alt':bg_alt,'datatype':datatype,
                'winsize':winsize,'wavelet':wavelet,'noisethresh':noisethresh,'depolnoisethresh':depolnoisethresh,
                'cloudthresh':cloudthresh,'CWTwidth':CWTwidth,'minwidth':minwidth,'layerCWTrange':layerCWTrange,
                'doPBL':doPBL,'PBLwavelet':PBLwavelet,'PBLCWTrange':PBLCWTrange,'PBLwidth':PBLwidth,
                'PBLnumprofs':PBLnumprofs,'PBLwinsize':PBLwinsize,'maxaeroalt':maxaeroalt,
                'sigma0':sigma0,'depolsigma0':depolsigma0,'waterthresh':waterthresh,'icethresh':icethresh,
                'smokethresh':smokethresh,'dustthresh':dustthresh,'refalt':refalt,'PBLrat':PBLrat,'molrat':molrat,'moldepol':moldepol,
                'method':method,'lrat':lrat,'mode':mode,'NRB_limits':NRB_limits,'depol_limits':depol_limits,
                'back_limits':back_limits,'ext_limits':ext_limits,'hours':hours,'fsize':fsize,'interpolate':interpolate}
                
    for f in rawfiles:
        newdate=MPLtodatetime(f)
        if newdate>=startdate and newdate<=enddate:
            if chunksize=='days':
                if olddate.date()==newdate.date():
                    tempfilelist.append(f)
                    if f==rawfiles[-1]:
                        if verbose:
                            print "Processing file group ending on {0}-{1}-{2}".format(olddate.year,olddate.month,olddate.day)
                            print len(tempfilelist)
                        prockwargs['rawfiles']=tempfilelist
                        mpltemp=fileproc(**prockwargs)
                        MPLlist.append(mpltemp)
                    olddate=newdate
                else:
                    if verbose:
                        print "Processing file group ending on {0}-{1}-{2}".format(olddate.year,olddate.month,olddate.day)
                        print len(tempfilelist)
                    prockwargs['rawfiles']=tempfilelist
                    mpltemp=fileproc(**prockwargs)
                    MPLlist.append(mpltemp)
                    tempfilelist=[f]
                    olddate=newdate
            elif chunksize=='weeks':
                if newdate-olddate<=datetime.timedelta(7):
                    tempfilelist.append(f)
                    if f==rawfiles[-1]:
                        prockwargs['rawfiles']=tempfilelist
                        mpltemp=fileproc(**prockwargs)
                        MPLlist.append(mpltemp)
                    olddate=newdate
                else:
                    if verbose:
                        print "Processing file group ending on {0}-{1}-{2}".format(olddate.year,olddate.month,olddate.day)
                    prockwargs['rawfiles']=tempfilelist
                    mpltemp=fileproc(**prockwargs)
                    MPLlist.append(mpltemp)
                    tempfilelist.append(f)
                    olddate=newdate
            elif chunksize=='months':
                if olddate.month==newdate.month:
                    tempfilelist.append(f)
                    if f==rawfiles[-1]:
                        prockwargs['rawfiles']=tempfilelist
                        mpltemp=fileproc(**prockwargs)
                        MPLlist.append(mpltemp)
                    olddate=newdate
                else:
                    if verbose:
                        print "Processing file group ending on {0}-{1}-{2}".format(olddate.year,olddate.month,olddate.day)                    
                    prockwargs['rawfiles']=tempfilelist
                    mpltemp=fileproc(**prockwargs)
                    MPLlist.append(mpltemp)
                    tempfilelist.append(f)
                    olddate=newdate
    os.chdir(olddir)
    return MPLlist
        
        
    
def quickplot(filename, datadir = None, savefigs = False):
    
#    MPLdat = mtools.MPL()
#    
#    MPLdat.fromHDF(filename[0])
    
    NRBcopol = pan.read_hdf(filename[0],'LNC')
    NRBdepol = pan.read_hdf(filename[0],'MPL')    
    
    if savefigs:    
        if datadir is None:
            datadir = mtools.set_dir("Select Folder for Figures to be depositied")
        
        os.chdir(datadir)
        
        try:
            os.chdir('Figures')
        except WindowsError:
            os.makedirs('Figures')
            os.chdir('Figures')
    
    #create figure and plot image of depolarization ratios
    fsize = 18 #baseline font size
    ar = 2.0  #aspect ratio
    figheight = 12 #inches
    
    plt.rc('font', family='serif', size=fsize)
    
       
    fig = plt.figure()
    fig.clf()
    h_set = range(1,25)
    h_set = map(str,h_set)
    
#    NRBcopol = MPLdat.NRB[0]
#    
#    NRBdepol = MPLdat.depolrat[0]
#    
    datetime = NRBcopol.index
    alt = NRBcopol.columns
    
    print 'Generating Figure'
    
    ax1 = fig.add_subplot(2,1,1)
    im1 = mplot.depol_plot(fig, ax1, ar, datetime,alt[::-1],NRBcopol.T[::-1], (0,.5), fsize = fsize)
    cbar1 = fig.colorbar(im1, orientation = 'vertical', aspect = 6)
    cbar1.set_ticks(np.arange(0,0.6,0.1))
    cbar1.set_ticklabels(np.arange(0,0.6,0.1))
    cbar1.ax.tick_params(labelsize = fsize)    
#    cbar1 = fig.colorbar(im1, orientation = 'vertical', aspect = 6)
#    cbar1.ax.tick_params(labelsize = fsize)
#    cbar1.ax.set_ylabel('$[km^{-1}sr^{-1}]$')
    mplot.dateticks(ax1, datetime, fsize = fsize, tcolor = 'w')
    ax1.set_xticklabels([])
    t1 = ax1.set_title('CORALNet Depol Ratio', fontsize = fsize+10)
    t1.set_y(1.03)
            
    ax2 = fig.add_subplot(2,1,2)
    im2 = mplot.depol_plot(fig, ax2, ar, datetime,alt[::-1],NRBdepol.T[::-1], (0,0.5), fsize = fsize)
    cbar2 = fig.colorbar(im2, orientation = 'vertical', aspect = 6)
    cbar2.set_ticks(np.arange(0,0.6,0.1))
    cbar2.set_ticklabels(np.arange(0,0.6,0.1))
    cbar2.ax.tick_params(labelsize = fsize)
    #set axis ranges and tickmarks based on data ranges
    mplot.dateticks(ax2, datetime, fsize = fsize)
    ax2.set_xlabel('Time [Local]',fontsize = fsize+4)
    fig.autofmt_xdate()
    t2 = ax2.set_title('miniMPL Depol Ratio',fontsize = fsize+10)
    t2.set_y(1.03)   
    
    plt.show()
    
    if savefigs:    
        ##plt.savefig(savetitle,dpi = 100, edgecolor = 'b', bbox_inches = 'tight')
        fig.set_size_inches(figheight*ar,figheight) 
        plt.savefig(filename+'NRB.png')
    print 'Done'
    
def MPLtodatetime(filename):
       
    year  = int(filename[0:4])
    month = int(filename[4:6])
    day   = int(filename[6:8])
    hour  = int(filename[8:10])
    
    dt_out = datetime.datetime(year,month,day,hour)
    return dt_out

def H5todatetime(filename):
    datetxt=filename.split('_proc.h5')
    [starttxt,endtxt]=datetxt[0].split('-')
    
    dt_start=MPLtodatetime(starttxt)
    if endtxt:
        dt_end=MPLtodatetime(endtxt)
    else:
        dt_end=[]
    
    return [dt_start,dt_end]
    

def datetimetoMPL(dt_in):
    
    year  = str(dt_in.year)
    month = str(dt_in.month).zfill(2)
    day   = str(dt_in.day).zfill(2)
    hour  = str(dt_in.hour).zfill(2)
    
    filename = '{0}{1}{2}{3}00.mpl'.format(year,month,day,hour)
    return filename

def datetimetoH5(dt_range):
    
    [dt_start,dt_end]=dt_range
    
    startyear  = str(dt_start.year)
    startmonth = str(dt_start.month).zfill(2)
    startday   = str(dt_start.day).zfill(2)
    starthour  = str(dt_start.hour).zfill(2)
    starttxt = '{0}{1}{2}{3}00'.format(startyear,startmonth,startday,starthour)
    
    if dt_end:
        endyear  = str(dt_end.year)
        endmonth = str(dt_end.month).zfill(2)
        endday   = str(dt_end.day).zfill(2)
        endhour  = str(dt_end.hour).zfill(2)
        endtxt = '{0}{1}{2}{3}00'.format(endyear,endmonth,endday,endhour)
        filename='{0}-{1}_proc.h5'.format(starttxt,endtxt)
    else:
        filename='{0}_proc.h5'.format(starttxt)
    
    return filename
        
def clearall():
    """clear all globals"""
    for uniquevar in [var for var in globals().copy() if var[0] != "_" and var != 'clearall']:
        del globals()[uniquevar]


        
if __name__ == '__main__':
    
    os.chdir(datalib)
#    os.chdir('C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Ucluelet Files')
#    olddir=os.getcwd()
#    datafile=mtools.set_dir('Select folder containing .mpl data to process')
#    os.chdir(datafile)
    altrange=np.arange(0.150,10.030,0.030)
    timestep='240S'
    savetype='standard'
    procsavepath=proclib
    plotsavepath=figurelib
    startdate=datetime.datetime(2011,05,03,00)
    enddate=datetime.datetime(2041,05,03,15)
    recalc=False
    
    NRBmask=False
    NRBthresh=3.0
    NRBmin=0.5
    NRBminalt=0.150
    NRBnumprofs=3
    SNRmask=True
    SNRthresh=0.5
    
    doPBL=False
    PBLCWTrange=np.arange(2,15)
    PBLwidth=2
    PBLnumprofs=3
    PBLwinsize=3
    
    
    PBLrat=30.0
    molrat=0.0
    moldepol=0.0035
    
    molthresh=1.0
    noisethresh=1.0
    depolnoisethresh=1.0
    sigma0=0.01
    depolsigma0=0.05
    cloudthresh=(1.0,0.20)
    waterthresh=0.10
    icethresh=0.35
    smokethresh=0.10
    dustthresh=0.20
    
    docorrection=True
    cor_numprofs=1
    cor_winsize=1
 
    doplot=True
    saveplot=True
    showplot=True
    dolayers=True
    
    dolayerplot=True
    docorplot=True
    verbose=False
#    hours=['02','04','06','08','10','12','14','16','18','20','22']
#    hours=['06','12','18']
    hours=['12']
    NRB_limits=(0.0,0.2,0.05) 
    depol_limits=(0.0,0.5,0.1)
    back_limits=(0.0,2e-4,5e-5)
    ext_limits=(0.0,1e-2,2e-3)
    interpolate='none'
    
    kwargs={'recalc':recalc,'altrange':altrange,'timestep':timestep,'savetype':savetype,'procsavepath':procsavepath,
            'plotsavepath':plotsavepath,'startdate':startdate,'enddate':enddate,
            'SNRthresh':SNRthresh,'molthresh':molthresh,'noisethresh':noisethresh,'depolnoisethresh':depolnoisethresh,
            'doplot':doplot,'saveplot':saveplot,'showplot':showplot,'dolayers':dolayers,
            'docorrection':docorrection,'cor_numprofs':cor_numprofs,'cor_winsize':cor_winsize,
            'dolayerplot':dolayerplot,'docorplot':docorplot,
            'verbose':verbose,'NRBmask':NRBmask,'SNRmask':SNRmask,'NRB_limits':NRB_limits,
            'depol_limits':depol_limits,'back_limits':back_limits,'ext_limits':ext_limits,
            'interpolate':interpolate, 'sigma0':sigma0, 'depolsigma0':depolsigma0,
            'PBLrat':PBLrat,'molrat':molrat,'moldepol':moldepol,
            'cloudthresh':cloudthresh,'waterthresh':waterthresh,'icethresh':icethresh,
            'smokethresh':smokethresh,'dustthresh':dustthresh,'hours':hours,
            'doPBL':doPBL,'PBLCWTrange':PBLCWTrange,'PBLwidth':PBLwidth,'PBLnumprofs':PBLnumprofs,
            'PBLwinsize':PBLwinsize,'newdir':datalib}

    mpl1=fileproc(**kwargs)
#    processall(**kwargs)
#    os.chdir(olddir)
    