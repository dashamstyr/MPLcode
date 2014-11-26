import numpy as np
import os,glob
import MPLtools as mtools
import MPLprocesstools as mproc
import datetime
import matplotlib.pyplot as plt
import MPLplot as mplot
import pandas as pan

def fileproc(newdir, **kwargs):

    """
    ----------------------------------------------------------------------------
    Uses tools created in MPL_tools to open all files in a folder and resample
    them to a regular spacing in altitude/date the concatenates them into one
    pandas dataframe and saves it as a 
    created: July 05, 2012
    updated: May 09,2014
    ----------------------------------------------------------------------------     
    inputs:
        newdir = directory containing .mpl files
    **kwargs:
         altrange = a list of altitudes to process(becomes pandas columns)
         starttime = time point for first profile
         endtime = time point for last profile
         timestep = time between profiles
         interactive = boolean determining if user will select files or if all in folder will be processed
         doplot = boolean determining if plot will be made from processed file
         saveplot = boolean determining if plot will be saved
         verbose = boolean determining if messages will be displayed
         
    
    """
    
    altrange = kwargs.get('altrange',np.arange(150,15000,30))        
    timestep = kwargs.get('timestep','60S')
    interactive = kwargs.get('interactive',True)
    doplot = kwargs.get('doplot',False)
    saveplot = kwargs.get('saveplot',False)
    savetype = kwargs.get('savetype','standard')
    showplot = kwargs.get('showplot',False)
    NRBmask = kwargs.get('NRBmask',False)
    SNRmask = kwargs.get('SNRmask',False)
    verbose = kwargs.get('verbose',False)
    NRB_limits = kwargs.get('NRB_limits',(0.0,1.0,0.2))  
    depol_limits = kwargs.get('depol_limits',(0.0,0.5,0.1))
    
    #starttime and endtime are defined later   
    
    olddir = os.getcwd()
    os.chdir(newdir)
    
    if interactive:
        rawfiles = mtools.get_files('Select MPL files', filetype = ('.mpl', '*.mpl'))
    else:
        rawfiles = glob.glob('*.mpl')
      
    #open, altitude resample, and concatenate data and mask files
    rawfiles.sort()
    
    if verbose:
        print rawfiles
    
    [path,startfile] = os.path.split(rawfiles[0])
    [path,endfile] = os.path.split(rawfiles[-1])
    
    starttime = kwargs.get('starttime',MPLtodatetime(startfile))
    endtime = kwargs.get('endtime',MPLtodatetime(endfile))
    
    if savetype=='standard':
        if len(rawfiles) == 1:   
            d_filename = '{0}_proc.h5'.format(startfile.split('.')[0])
        else:        
            d_filename = '{0}-{1}_proc.h5'.format(startfile.split('.')[0],endfile.split('.')[0])
    elif savetype=='IDL':
        if len(rawfiles) == 1:   
            d_filename = '{0}_IDL.h5'.format(startfile.split('.')[0])
        else:        
            d_filename = '{0}-{1}_IDL.h5'.format(startfile.split('.')[0],endfile.split('.')[0])
    
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
    
    if NRBmask:
        MPLdat_event=mproc.NRB_mask_all(MPLdat_event)
    
    if SNRmask:
        MPLdat_event=mproc.SNR_mask_depol(MPLdat_event)
    
    if os.path.isdir('Processed'):
        os.chdir('Processed')
    else:
        os.makedirs('Processed')
        os.chdir('Processed')
    
    if verbose:
        print 'Saving '+d_filename
    
    if savetype=='standard':
        MPLdat_event.save_to_HDF(d_filename)
    elif savetype=='IDL':
        MPLdat_event.save_to_IDL(d_filename)
    
    if verbose:
        print 'Done'
    
    if doplot:
        import MPLplot as mplot 
        savefilename = '{0}.png'.format(d_filename.split('.')[0])
        savepath = os.path.join(newdir,'Figures')
        mplot.doubleplot(MPLdat_event,saveplot=saveplot,showplot=showplot,
                         verbose=verbose,savefilepath=savepath,savefilename=savefilename,
                         NRB_limits=NRB_limits,depol_limits=depol_limits)
    
    os.chdir(olddir)
    
    
def quickplot(filename, datadir = [], savefigs = False):
    
#    MPLdat = mtools.MPL()
#    
#    MPLdat.fromHDF(filename[0])
    
    NRBcopol = pan.read_hdf(filename[0],'LNC')
    NRBdepol = pan.read_hdf(filename[0],'MPL')    
    
    if savefigs:    
        if not datadir:
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
    
    import MPLtools as mtools
    
    altrange=range(150,15180,30)
    timestep='60S'
    os.chdir('C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Ucluelet Files')
    procdir = mtools.set_dir('Select folder to process files from')
    fileproc(procdir,doplot=True,saveplot=True,showplot=True,verbose=True,
             altrange=altrange,timestep=timestep,NRB_limits=(0.0,1.0,0.2),savetype='standard')
    
#    plotfile = mtools.get_files('Select File to Plot', filetype = ('.h5', '*.h5'))
#    quickplot(plotfile)
    