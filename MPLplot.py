import numpy as np
import os, sys
import numpy as np
import datetime as dt
import bisect
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pan
import MPLtools as mtools
import MPLprocesstools as mproc
from datetime import datetime

def custom_cmap(maptype,numvals,overcolor,undercolor):
    if maptype=="customjet":
        cdict = {'red': ((0.00, 0, 0),
                         (0.35, 0, 0),
                         (0.66, 1, 1),
                         (0.89,1, 1),
                         (1, 0.5, 0.5)),
             'green': ((0.00, 0, 0),
                       (0.125,0, 0),
                       (0.375,1, 1),
                       (0.64,1, 1),
                       (0.91,0,0),
                       (1, 0, 0)),
             'blue': ((0.00,0.5,0.5),
                      (0.11, 1, 1),
                      (0.34, 1, 1),
                      (0.65,0, 0),
                      (1, 0, 0))}        
    elif maptype=="greyscale":
        cdict = {'red':   ((0.0,0,0),
                           (1.0,1,1)),
                 'green': ((0.0,0,0),
                           (1.0,1,1)),
                 'blue':  ((0.0,0,0),
                           (1.0,1,1))}

    my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,numvals)
    my_cmap.set_over(overcolor)
    my_cmap.set_under(undercolor)        
    
    return my_cmap
    
def depol_plot(fig,ax,data,xdata,ydata,**kwargs):    
    #set colormap to be the same as 'jet' with the addition of white color for
    #depol ratios set to identically zero because they couldn't be calculated
    ar=kwargs.get('ar',2.0)
    vrange=kwargs.get('vrange',[0,1])
    fsize=kwargs.get('fsize',21)
    maptype=kwargs.get('maptype','customjet')
    orientation=kwargs.get('orientation','Vertical')
    overcolor=kwargs.get('overcolor','w')
    undercolor=kwargs.get('undercolor','k')
    numvals=kwargs.get('numvals',50)
    
    my_cmap=custom_cmap(maptype=maptype,numvals=numvals,overcolor=overcolor,undercolor=undercolor)    
    im = ax.imshow(data, vmin=vrange[0], vmax=vrange[1], cmap = my_cmap)
    forceAspect(ax,ar)        
    altticks(ax, ydata, fsize=fsize)  
    if orientation=='vertical':
        ax.set_ylabel('Altitude [m]', fontsize = fsize+4, labelpad = 15)
    elif orientation=='horizontal':
        ax.set_ylabel('Horizontal Range [m]', fontsize = fsize+4, labelpad = 15)
    for line in ax.yaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(1)        
    ax.axis('tight')

    return im

def backscatter_plot(fig, ax, data, xdata, ydata, **kwargs):    
    #set colormap to be the same as 'jet' with the addition of white color for
    #depol ratios set to identiacally zero because they couldn't be calculated
    ar=kwargs.get('ar',2.0)
    vrange=kwargs.get('vrange',[0,1])
    fsize=kwargs.get('fsize',21)
    maptype=kwargs.get('maptype','customjet')
    orientation=kwargs.get('orientation','Vertical')
    overcolor=kwargs.get('overcolor','w')
    undercolor=kwargs.get('undercolor','k')
    numvals=kwargs.get('numvals',50)
    
    my_cmap=custom_cmap(maptype=maptype,numvals=numvals,overcolor=overcolor,undercolor=undercolor)  
    im = ax.imshow(data, vmin=vrange[0], vmax=vrange[1], cmap = my_cmap)
    forceAspect(ax,ar)       
    altticks(ax, ydata, fsize = fsize, tcolor = 'w')
    
    if orientation=='vertical':
        ax.set_ylabel('Altitude [m]', fontsize = fsize+4, labelpad = 15)
    elif orientation=='vorizontal':
        ax.set_ylabel('Horizontal Range [m]', fontsize = fsize+4, labelpad = 15)
    for line in ax.yaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(1)        
    ax.axis('tight')

    return im
    

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def dateticks(ax, axisdat,hours = [], fsize = 21, tcolor = 'k'):
    
    dold = axisdat[0].strftime('%d')
    hold = axisdat[0].strftime('%H')
    tickmarks = []
    ticklabels = []
    n = 0

    
    for d in axisdat:
        dtemp = d.strftime('%d')
        if dtemp != dold:
            ticklabels.append(d.strftime('%b %d'))
            tickmarks.append(n)
        else:
            htemp = d.strftime('%H')
            if not hours:
                if htemp != hold:
                    ticklabels.append(d.strftime('%H'))
                    tickmarks.append(n)
            else:
                if htemp in hours and htemp != hold:
                    ticklabels.append(d.strftime('%H'))
                    tickmarks.append(n)
            hold = htemp

        dold = dtemp
        n += 1
    
    ax.set_xticks(tickmarks)
    ax.set_xticklabels(ticklabels, fontsize = fsize)

    for line in ax.xaxis.get_ticklines():
        line.set_color(tcolor)
        line.set_markersize(10)
        line.set_markeredgewidth(2)

def altticks(ax, axisdat, numticks = 5, fsize = 21, tcolor = 'k'):

    numpoints = len(axisdat)
    step = numpoints//numticks
    tickmarks = range(0,numpoints,step)
    ticklabels = [str(int(t)) for t in axisdat[::step]]

    ax.set_yticks(tickmarks)
    ax.set_yticklabels(ticklabels, fontsize = fsize)

    for line in ax.yaxis.get_ticklines():
        line.set_color(tcolor)
        line.set_markersize(10)
        line.set_markeredgewidth(3)

def vertprof(df, altrange, exact_times, plot_type = 'line', zeromask = False, 
             savefig = False, filename = 'tempfig.png'):    
    
    minalt = altrange[0]
    maxalt = altrange[1]
    
    daterange = df.index
    ymin = df.columns[0]
    ymax = df.columns[-1]
    
    if ymax > maxalt:
        df = df.loc[:,:maxalt]
    
    if minalt > ymin:
        df.loc[:,:minalt] = 'nan'
    
    approx_times = []
    
   
    for ts in exact_times:    
        i = bisect.bisect_left(daterange, ts)
        approx_times.append(min(daterange[max(0, i-1): i+2], key=lambda t: abs(ts - t)))
    
    plt.clf()
    fig = plt.figure()
    
    numprof = len(approx_times)
    
    for n in range(numprof): 
        print approx_times[n]
        s = df.ix[approx_times[n]]
        
        if zeromask:  
            s = s[s>0]
        else:
            zeroline = np.zeros_like(s)
            
        alt = s.index
        print s.max()
        ax = fig.add_subplot(1,numprof,n+1)
        
        if plot_type == 'line':
            im = ax.plot(s,alt, linewidth = 4)
            plt.ylim([ymin,ymax])
        
        if plot_type == 'scatter':
            im = ax.scatter(s,alt)
            plt.ylim([ymin,ymax])
        
        try:
            zeroline 
        except NameError:
            continue
        else:
            ax.plot(zeroline,alt,'r--', linewidth = 2)
            
        plt.yticks(fontsize = 21)
        plt.ylabel('Altitude [m]', fontsize = 21)
        
        for line in ax.yaxis.get_ticklines():
            line.set_color('k')
            line.set_markersize(6)
            line.set_markeredgewidth(2)
            
        for line in ax.xaxis.get_ticklines():
            line.set_color('k')
            line.set_markersize(6)
            line.set_markeredgewidth(2)    
            
        plt.xticks(fontsize = 21)
        plt.title(approx_times[n], fontsize = 21)
        fig.subplots_adjust(wspace = 0.5)
    
    plt.show()
    
    if savefig:
        plt.savefig(filename)        
    
def doubleplot(datafile,**kwargs):    
    """
    inputs
    --------------------------------------------------------------------------
    datafile = processed .h5 file or MPL class object to plot

    kwargs:
    altrange = range of altitude values to plot
    timestep = string represtning step between profiles
    starttime = datetime object denoting time of first profile
    endtime = datetime object denoting time of last profile
    hours = list of strings denoting times to label in plot
    fsize = baseline font size
    ar = aspect ratio
    figheight = figure height in inches
    NRB_limits = (min,max,step) defines color scale for NRB
    depol_limits = (min,max,step) defines color scale for depol
    dpi = plot resolutoion in dots per inch
    savefile = boolean to determine whather plot will be saved
    savefilepath = location to save file
    showplot = boolean to determine whether to display plot
    verbose = boolean to determine whether to display messages
    SNRmask = boolean to determine whether to apple SNR masking
    SNRthresh = SNR threshold to apply (defaults to 1)
    
    """
    
    #set kwarg defaults 
    starttime = kwargs.get('starttime',[])
    endtime = kwargs.get('endtime',[])
    timestep = kwargs.get('timestep',[])
    altrange = kwargs.get('altrange',[])
    hours = kwargs.get('hours',['03','06', '09','12', '15','18','21']) 
    fsize = kwargs.get('fsize',18) #baseline font size
    ar = kwargs.get('ar',2.0)  #aspect ratio
    figheight = kwargs.get('figheight',12) #inches    
    NRB_limits = kwargs.get('NRB_limits',(0.0,1.0,0.2))  
    depol_limits = kwargs.get('depol_limits',(0.0,0.5,0.1))
    dpi = kwargs.get('dpi',100)
    saveplot = kwargs.get('saveplot',True)
    colormap=kwargs.get('colormap','customjet')
    orientation=kwargs.get('orientation','vertical')
    
    if type(datafile)==str:
        savefilename = kwargs.get('savefilename','{0}.png'.format(datafile.split('.')[0]))
    else:
        savefilename = kwargs.get('savefilename','MPLdoubleplot.png')
    savefilepath = kwargs.get('savefilepath',os.path.join(os.getcwd(),'Figures'))
    showplot = kwargs.get('showplot',False)
    verbose = kwargs.get('verbose',False)
    SNRmask = kwargs.get('SNRmask',False)
    SNRthresh = kwargs.get('SNRthresh',3) 
    
    NRB_min = NRB_limits[0]
    NRB_max = NRB_limits[1]
    NRB_step = NRB_limits[2]
    
    depol_min = depol_limits[0]
    depol_max = depol_limits[1]
    depol_step = depol_limits[2]
    
    if type(datafile) in [str,unicode]:
        MPLevent = mtools.MPL()
        
        MPLevent.fromHDF(datafile)    
    else:
        MPLevent = datafile

    if not MPLevent.depolrat:
        MPLevent.calculate_depolrat()
        
    if len(altrange)>0:
        MPLevent.alt_resample(altrange)    
    
    if [i for i in [starttime,endtime,timestep] if len(i)>0]:
        MPLevent.time_resample(timestep=timestep,starttime=starttime,endtime=endtime)
        
      
    #Now generate a new figure.
    
    if SNRmask:
        MPLevent=mproc.SNR_mask_depol(MPLevent)
    
    copol = MPLevent.NRB[0]
    depol = MPLevent.depolrat[0]
        
    copol.fillna(-99999,inplace=True)
    depol.fillna(-99999,inplace=True)
    alt = copol.columns
    datetime = copol.index    
    
    #create figure and plot image of depolarization ratios    
    plt.rc('font', family='serif', size=fsize)
    
    
    fig = plt.figure(0)
    
    h_set = range(1,25)
    h_set = map(str,h_set)
    
    if verbose:
        print 'Generating Figure'
        
    ax1 = fig.add_subplot(2,1,1)
    im1 = backscatter_plot(fig, ax1,copol.T[::-1],datetime,alt[::-1],ar=ar, 
                           vrange=(NRB_min,NRB_max),fsize=fsize,maptype=colormap,
                            orientation=orientation)
    cbar1 = fig.colorbar(im1, orientation = 'vertical', aspect = 6, extend='both')
    cbar1.set_ticks(np.arange(NRB_min,NRB_max+NRB_step,NRB_step))
    cbar1.set_ticklabels(np.arange(NRB_min,NRB_max+NRB_step,NRB_step))
    cbar1.ax.tick_params(labelsize = fsize)
    cbar1.ax.set_ylabel('$[counts*km^{2}/(\mu s*\mu J)$')
    dateticks(ax1, datetime, hours = hours, fsize = fsize, tcolor = 'w')
    ax1.set_xticklabels([])
    t1 = ax1.set_title('Normalized Relative Backscatter', fontsize = fsize+10)
    t1.set_y(1.03)
            
    ax2 = fig.add_subplot(2,1,2)
    im2 = depol_plot(fig, ax2, depol.T[::-1],datetime,alt[::-1],ar=ar, 
                     vrange=(depol_min,depol_max),fsize=fsize,maptype=colormap,
                        orientation=orientation)
    cbar2 = fig.colorbar(im2, orientation = 'vertical', aspect = 6, extend='both')
    cbar2.set_ticks(np.arange(depol_min,depol_max+depol_step,depol_step))
    cbar2.set_ticklabels(np.arange(depol_min,depol_max+depol_step,depol_step))
    cbar2.ax.tick_params(labelsize = fsize)
    #set axis ranges and tickmarks based on data ranges
    dateticks(ax2, datetime, hours = hours, fsize = fsize)
    ax2.set_xlabel('Time [Local]',fontsize = fsize+4)
    fig.autofmt_xdate()
    t2 = ax2.set_title('Linear Depolarization Ratio',fontsize = fsize+10)
    t2.set_y(1.03)    
    
    
    fig.set_size_inches(figheight*ar,figheight) 
    if saveplot:
        olddir=os.getcwd()
        if os.path.isdir(savefilepath):
            os.chdir(savefilepath)
        else:
            os.mkdir(savefilepath)
            os.chdir(savefilepath)
        fig.canvas.print_figure(savefilename,dpi = dpi, edgecolor = 'b', bbox_inches = 'tight') 
        os.chdir(olddir)
    if showplot:
        fig.canvas.draw()
    del fig, cbar1,cbar2,MPLevent
    
    if verbose:
        print 'Done'
#    del MPLevent
    
if __name__=='__main__':       
    altrange = np.arange(150,15030,30)
    starttime = []
    endtime = []
    timestep = '120S'
    hours = ['03','06','09','12','15','18','21']
    depol_limits=(0.0,0.5,0.1)
    NRB_limits=(0.0,1.0,0.2)

    os.chdir('C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA\Ucluelet Files\Processed')

    filepath = mtools.get_files('Select Processed MPL file(s)', filetype = ('.h5', '*.h5'))
    
    if len(filepath)==1:
        [d_path,d_filename] = os.path.split(filepath[0])    
        savefilename = '{0}.png'.format(d_filename.split('.')[0])
    else:
        [d_path,startfile] = os.path.split(filepath[0])
        [d_path,endfile] = os.path.split(filepath[-1])
        savefilename = '{0}-{1}.png'.format(startfile.split('.')[0],endfile.split('.')[0])          
    savepath = os.path.join(os.path.split(d_path)[0],'Figures')

    for f in filepath:  
        MPLtemp = mtools.MPL()
        MPLtemp.fromHDF(f)        
        MPLtemp.alt_resample(altrange)
        
        try:
            MPLevent.append(MPLtemp)
        except NameError:
            MPLevent = MPLtemp
   
    MPLevent.header.sort_index(inplace=True)
    
    for n in range(MPLevent.header['numchans'][0]):
        MPLevent.data[n].sort_index(inplace=True)
    
    MPLevent.time_resample(timestep,starttime,endtime)  
    MPLevent.calc_all()
    MPLevent=mproc.NRB_mask_all(MPLevent) 
    
    kwargs = {'saveplot':True,'showplot':True,'verbose':True,
                'savefilepath':savepath,'savefilename':savefilename,
                'hours':hours,'depol_limits':depol_limits,'NRB_limits':NRB_limits,'SNRmask':True,
                'colormap':'customjet','orientation':'vertical'}
    
    doubleplot(MPLevent,**kwargs)
    
    del MPLevent
    
