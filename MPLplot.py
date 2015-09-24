#from __future__ import absolute_import
import os,sys,site

import numpy as np
import os, sys
import numpy as np
import bisect
import operator
from copy import deepcopy
import datetime
import MPLtools as mtools

if sys.platform=='win32':
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors    
    import matplotlib.ticker as mtick
    from mpl_toolkits.axes_grid1 import make_axes_locatable
else:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    from matplotlib import colors as colors
    from matplotlib import ticker as mtick
    from mpl_toolkits.axes_grid1 import make_axes_locatable

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
    
def top_plot(ax, data, xdata, ydata, **kwargs):    

    ar=kwargs.get('ar',2.0)
    vrange=kwargs.get('vrange',[0,1])
    fsize=kwargs.get('fsize',21)
    maptype=kwargs.get('maptype','customjet')
    orientation=kwargs.get('orientation','Vertical')
    overcolor=kwargs.get('overcolor','w')
    undercolor=kwargs.get('undercolor','k')
    numvals=kwargs.get('numvals',50)
    interpolation=kwargs.get('interpolation',None)
    
    my_cmap=custom_cmap(maptype=maptype,numvals=numvals,overcolor=overcolor,undercolor=undercolor)  
    im = ax.imshow(data, vmin=vrange[0], vmax=vrange[1], cmap = my_cmap, interpolation=interpolation)
    forceAspect(ax,ar)       
    altticks(ax, ydata, fsize = fsize, tcolor = 'w')
    
    if orientation=='vertical':
        ax.set_ylabel('Altitude [km]', fontsize = fsize+4, labelpad = 15)
    elif orientation=='vorizontal':
        ax.set_ylabel('Horizontal Range [km]', fontsize = fsize+4, labelpad = 15)
    for line in ax.yaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(1)        
    ax.axis('tight')

    return im

def bottom_plot(fig,ax,data,xdata,ydata,**kwargs):    

    ar=kwargs.get('ar',2.0)
    vrange=kwargs.get('vrange',[0,1])
    fsize=kwargs.get('fsize',21)
    maptype=kwargs.get('maptype','customjet')
    orientation=kwargs.get('orientation','Vertical')
    overcolor=kwargs.get('overcolor','w')
    undercolor=kwargs.get('undercolor','k')
    numvals=kwargs.get('numvals',50)
    interpolation=kwargs.get('interpolation',None)
    my_cmap=custom_cmap(maptype=maptype,numvals=numvals,overcolor=overcolor,undercolor=undercolor)    
    im = ax.imshow(data, vmin=vrange[0], vmax=vrange[1], cmap = my_cmap, interpolation=interpolation)
    forceAspect(ax,ar)        
    altticks(ax, ydata, fsize=fsize)  
    if orientation=='vertical':
        ax.set_ylabel('Altitude [km]', fontsize = fsize+4, labelpad = 15)
    elif orientation=='horizontal':
        ax.set_ylabel('Horizontal Range [km]', fontsize = fsize+4, labelpad = 15)
    for line in ax.yaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(1)        
    ax.axis('tight')

    return im    

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def dateticks(ax, axisdat,hours = None, fsize = 21, tcolor = 'k'):
    
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
            if hours is None:
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
    ticklabels = ["%.2F"%(t) for t in axisdat[::step]]

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
    
    numfigs=len(plt.get_fignums())
    fig = plt.figure(numfigs+1)
    
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
        plt.ylabel('Altitude [km]', fontsize = 21)
        
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

def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)
    
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
    topplot_limits = (min,max,step) defines color scale for NRB
    bottomplot_limits = (min,max,step) defines color scale for depol
    dpi = plot resolutoion in dots per inch
    savefile = boolean to determine whather plot will be saved
    savefilepath = location to save file
    showplot = boolean to determine whether to display plot
    verbose = boolean to determine whether to display messages
    SNRmask = boolean to determine whether to apple SNR masking
    SNRthresh = SNR threshold to apply (defaults to 1)
    
    """
    
    #set kwarg defaults 
    starttime = kwargs.get('starttime',None)
    endtime = kwargs.get('endtime',None)
    timestep = kwargs.get('timestep',None)
    altrange = kwargs.get('altrange',None)
    hours = kwargs.get('hours',['03','06', '09','12', '15','18','21']) 
    fsize = kwargs.get('fsize',18) #baseline font size
    ar = kwargs.get('ar',2.0)  #aspect ratio
    figheight = kwargs.get('figheight',12) #inches    
    topplot_limits = kwargs.get('topplot_limits',(0.0,1.0,0.2))  
    bottomplot_limits = kwargs.get('bottomplot_limits',(0.0,0.5,0.1))
    dpi = kwargs.get('dpi',100)
    saveplot = kwargs.get('saveplot',True)
    colormap=kwargs.get('colormap','customjet')
    orientation=kwargs.get('orientation','vertical')
    toptype=kwargs.get('toptype','NRB')
    bottomtype=kwargs.get('bottomtype','depol')
    interpolation=kwargs.get('interpolation',None)
    
    if type(datafile)in [str,unicode]:
        savefilename = kwargs.get('savefilename','{0}.png'.format(datafile.split('.')[0]))
    else:
        savefilename = kwargs.get('savefilename','MPLdoubleplot.png')
    savefilepath = kwargs.get('savefilepath',os.path.join(os.getcwd(),'Figures'))
    showplot = kwargs.get('showplot',False)
    verbose = kwargs.get('verbose',False)
    SNRmask = kwargs.get('SNRmask',False)
    SNRthresh = kwargs.get('SNRthresh',3)
    SNRtype = kwargs.get('SNRtype','NRB')
    
    topplot_min = topplot_limits[0]
    topplot_max = topplot_limits[1]
    topplot_step = topplot_limits[2]
    
    bottomplot_min = bottomplot_limits[0]
    bottomplot_max = bottomplot_limits[1]
    bottomplot_step = bottomplot_limits[2]
    
    if type(datafile) in [str,unicode]:
        MPLevent = mtools.MPL()
        
        MPLevent.fromHDF(datafile)    
    else:
        MPLevent = deepcopy(datafile)
        
    if altrange is not None:
        MPLevent.alt_resample(altrange)    
    
    if any([starttime,endtime,timestep]):
        MPLevent.time_resample(timestep=timestep,starttime=starttime,endtime=endtime)
        
      
    #Now generate a new figure.
        
    if toptype=='NRB':
        topdat = MPLevent.NRB[0]
        toptitle='Normalized Relative Backscatter'
        topunits='${counts*km^{2}}/{\mu s*\mu J}$'
    elif toptype=='backscatter':
        topdat = MPLevent.backscatter[0]
        toptitle='Backscatter Coefficient'
        topunits='$km^{-1}sr^{-1}$'
    
    if bottomtype=='depol':
        if not MPLevent.depolrat:
            MPLevent.calculate_depolrat()
        if SNRmask:
            MPLevent_masked=mtools.SNR_mask_depol(MPLevent,SNRthreshold=SNRthresh)
            bottomdat = MPLevent_masked.depolrat[0]
        else:
            bottomdat = MPLevent.depolrat[0]
        bottomtitle='Linear Depolarization Ratio'
        bottomunits=''
        
        
    elif bottomtype=='extinction':
        bottomdat = MPLevent.extinction[0]
        bottomtitle='Extinction Coefficient'
        bottomunits='$km^{-1}$'
        
    topdat.fillna(-99999,inplace=True)
    bottomdat.fillna(-99999,inplace=True)
    alt = topdat.columns
    dtindex = topdat.index    
    
    #create figure and plot image of depolarization ratios    
    plt.rc('font', family='serif', size=fsize)
    
    if showplot:
        numfigs=len(plt.get_fignums())
        
        fig = plt.figure(numfigs+1)
    else:
        fig=plt.figure()
    
    h_set = range(1,25)
    h_set = map(str,h_set)
    
    if verbose:
        print 'Generating Figure'
        
    ax1 = fig.add_subplot(2,1,1)
    im1 = top_plot(ax1,topdat.T[::-1],dtindex,alt[::-1],ar=ar, 
                           vrange=(topplot_min,topplot_max),fsize=fsize,maptype=colormap,
                            orientation=orientation, interpolation=interpolation)
    cbar1 = fig.colorbar(im1, orientation = 'vertical', aspect = 6, extend='both')
    cbar1.set_ticks(np.arange(topplot_min,topplot_max+topplot_step,topplot_step))
    cbar1.set_ticklabels(np.arange(topplot_min,topplot_max+topplot_step,topplot_step))
    cbar1.ax.tick_params(labelsize = fsize)
    cbar1.ax.set_ylabel(topunits)
    dateticks(ax1, dtindex, hours = hours, fsize = fsize, tcolor = 'w')
    ax1.set_xticklabels([])
    t1 = ax1.set_title(toptitle, fontsize = fsize+10)
    t1.set_y(1.03)
            
    ax2 = fig.add_subplot(2,1,2)
    im2 = bottom_plot(fig, ax2, bottomdat.T[::-1],dtindex,alt[::-1],ar=ar, 
                     vrange=(bottomplot_min,bottomplot_max),fsize=fsize,maptype=colormap,
                        orientation=orientation, interpolation=interpolation)
    cbar2 = fig.colorbar(im2, orientation = 'vertical', aspect = 6, extend='both')
    cbar2.set_ticks(np.arange(bottomplot_min,bottomplot_max+bottomplot_step,bottomplot_step))
    cbar2.set_ticklabels(np.arange(bottomplot_min,bottomplot_max+bottomplot_step,bottomplot_step))
    cbar2.ax.tick_params(labelsize = fsize)
    cbar2.ax.set_ylabel(bottomunits)
    #set axis ranges and tickmarks based on data ranges
    dateticks(ax2, dtindex, hours = hours, fsize = fsize)
    ax2.set_xlabel('Time [Local]',fontsize = fsize+4)
    fig.autofmt_xdate()
    t2 = ax2.set_title(bottomtitle,fontsize = fsize+10)
    t2.set_y(1.03)    
    
    
    fig.set_size_inches(figheight*ar,figheight) 
    if saveplot:
        olddir=os.getcwd()
        if os.path.isdir(savefilepath):
            os.chdir(savefilepath)
        else:
            os.mkdir(savefilepath)
            os.chdir(savefilepath)
        fig.savefig(savefilename,dpi = dpi, edgecolor = 'b', bbox_inches = 'tight') 
        os.chdir(olddir)
    if showplot:
        fig.canvas.draw()
    del fig, cbar1,cbar2,MPLevent
    
    if verbose:
        print 'Done'

def doubleprof(prof1,prof2,rangecor=True,deltaplot=True):
    plotprofs=[]
    
    for p in [prof1,prof2]:
        if rangecor:
            alts=p.index()
            plotprofs.append(p*alts**2)
        else:
            plotprofs.append(p)
    
    if deltaplot:
        normprofs=[]
        for p in plotprofs:
            tempmean=p.mean()
            normprofs.append(p/tempmean())
        deltaprof=(normprofs[1]-normprofs[2])*100.0/normprofs[1] 
    
    numfigs=len(plt.get_fignums())
    fig=plt.figure(numfigs+1)
    ax1=fig.add_subplot(211)
    ax1a=plotprofs[0].plot()
    ax1b=plotprofs[1].plot(secondary_y=True)
    align_yaxis(ax1a,0,ax1b,0)
    ax2=fig.add_subplot(212)
    deltaprof.plot()

def colormask_plot(mplin,**kwargs):
    #set color codes for different layers
    hours=kwargs.get('hours',['00','06','12','18'])
    fontsize=kwargs.get('fontsize',24)
    altrange=kwargs.get('altrange',None)
    datetimerange=kwargs.get('datetimerange',None)
    SNRmask=kwargs.get('SNRmask',False)
    SNRthresh=kwargs.get('SNRthresh',3.0)
    SNRtype=kwargs.get('SNRtype','NRB')
    saveplot=kwargs.get('saveplot',False)
    showplot=kwargs.get('showplot',True)
    plotfilepath=kwargs.get('plotfilepath',None)
    plotfilename=kwargs.get('plotfilename','testmaskfig.png')
    dpi = kwargs.get('dpi',100)    
    colordict=kwargs.get('colordict',{'Clear Air':0,
                                       'PBL':1,
                                       'Ice Cloud':2,
                                       'Water Cloud':3,
                                       'Mixed Cloud':4,
                                       'Dust':5,
                                       'Polluted Dust':6,
                                       'Smoke/Urban':7,
                                       'Unidentified Aerosol':8,
                                       'Insufficient Signal':9})
    
    Clear=(176.0/255.0,244.0/255.0,230.0/255.0)
    PBL=(255.0/255.0,69.0/255.0,0.0/255.0)
    Ice=(255.0/255.0,255.0/255.0,255.0/255.0)
    Water=(0.0/255.0,0.0/255.0,205.0/255.0)
    Mixed=(186.0/255.0,85.0/255.0,211.0/255.0)
    Dust=(255.0/255.0,193.0/255.0,37.0/255.0)
    Dust_Smoke=(199.0/255.0,97.0/255.0,20.0/255.0)
    Smoke_Urban=(75.0/255.0,75.0/255.0,100.0/255.0)
    Unidentified=(192.0/255.0,192.0/255.0,192.0/255.0)
    Insufficient=(0.0/255.0,0.0/255.0,0.0/255.0)
                                       
    
    maskmap=colors.ListedColormap([Clear,PBL,Ice,Water,Mixed,Dust,Dust_Smoke,Smoke_Urban,Unidentified,Insufficient])
    if SNRmask:
        mplmasked=mtools.SNR_mask_scene(mplin,SNRthresh=SNRthresh,datatype=SNRtype,inplace=False)
        maskin=mplmasked.scenepanel[0]['colormask']
    else:
        maskin=deepcopy(mplin.scenepanel[0]['colormask'])
        
    if altrange is not None:
        maskin=maskin.loc[:,(maskin.columns>altrange[0]) & (maskin.columns<altrange[-1])]
    
    if datetimerange is not None:
        maskin=maskin[(maskin.index>datetimerange[0]) & (maskin.index<datetimerange[-1])]
    
    plotmask=maskin.T[::-1]+0.5
    plotmask=plotmask.astype('float')
    times=maskin.index
    alts=maskin.columns
    
    if showplot:
        numfigs=len(plt.get_fignums())
        fig=plt.figure(numfigs+1,figsize=(30,5),dpi=100)
    else:
        fig=plt.figure(figsize=(30,5))
    
#    ax1=plt.subplot2grid((1,9),(0,0),rowspan=1,colspan=7)
#    cax1=plt.subplot2grid((1,9),(0,7),rowspan=1,colspan=1) 
    ax1=fig.add_subplot(111)
    image=ax1.imshow(plotmask,cmap=maskmap,interpolation='none',vmin=0,vmax=10,aspect='auto')
#    forceAspect(ax1,2.0)
    divider=make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.15)
    
    dateticks(ax1, times, hours = hours,fsize=fontsize)
    ax1.set_xlabel('Hours [Local]',fontsize=fontsize+4)
    ax1.set_ylabel('Altitude [km]', fontsize=fontsize+4)
    altticks(ax1, alts[::-1], fsize = fontsize, tcolor = 'k')
    cbar1=fig.colorbar(image,cax=cax1,orientation='vertical')
    cbar_ticklocs=np.arange(len(colordict))+0.5
    cbar1.set_ticks(cbar_ticklocs)
    cbar1.ax.tick_params(bottom='off',top='off',labelsize=fontsize-4)
    
    sortedlabels=[s[0] for s in sorted(colordict.iteritems(), key=operator.itemgetter(1))]
    cbarlabels=cbar1.set_ticklabels(sortedlabels)
    plt.tight_layout()
    if saveplot:
        if plotfilepath is not None:
            if os.path.isdir(plotfilepath):
                savename=os.path.join(plotfilepath,plotfilename)
            else:
                os.mkdir(plotfilepath)
                savename=os.path.join(plotfilepath,plotfilename)
        else:
            savename=os.path.join(os.getcwd(),plotfilename)
        fig.savefig(savename,dpi = dpi, edgecolor = 'b', bbox_inches = 'tight') 
    
    if showplot:
        fig.canvas.draw()

def plotall(mplfile,**kwargs):
    time_resample=kwargs.get('time_resample',False)
    alt_resample=kwargs.get('alt_resample',False)
    altrange = kwargs.get('altrange',np.arange(150,15000,30))        
    timestep = kwargs.get('timestep','60S')
    starttime = kwargs.get('starttime',None)
    endtime = kwargs.get('endtime',None)
    doplot=kwargs.get('doplot',True)
    docorplot=kwargs.get('docorplot',True)
    dolayerplot=kwargs.get('dolayerplot',True)
    saveplots=kwargs.get('saveplots',True)
    showplots=kwargs.get('showplots',True)
    plotsavepath = kwargs.get('plotsavepath','..\Figures')
    NRBmask=kwargs.get('NRBmask',True)
    NRBmasktype = kwargs.get('NRBmasktype','profile')
    NRBthresh=kwargs.get('NRBthresh',3.0)
    NRBmin=kwargs.get('NRBmin',0.5)
    NRBminalt=kwargs.get('NRBminalt',0.150)
    NRBnumprofs=kwargs.get('NRBnumprofs',1)
    NRBwinsize=kwargs.get('NRBwinsize',3)
    NRBinplace=kwargs.get('NRBinplace',False)
    NRB_limits = kwargs.get('NRB_limits',(0.0,1.0,0.2))  
    depol_limits = kwargs.get('depol_limits',(0.0,0.5,0.1))
    back_limits = kwargs.get('back_limits',(0.0,1e-7,2e-8))
    ext_limits = kwargs.get('ext_limits',(0.0,2e-6,4e-7))
    hours = kwargs.get('hours',['06','12','18']) 
    fsize = kwargs.get('fsize',18) #baseline font size
    SNRmask = kwargs.get('SNRmask',False)
    SNRthresh=kwargs.get('SNRthresh',1.0)
    SNRtype=kwargs.get('SNRtype','NRB')
    interpolate=kwargs.get('interpolate','None') #other methods include none, bilinear, gaussian and hermite
    #starttime and endtime are defined later  
    
    MPLdat_event=mtools.MPL()
    MPLdat_event.fromHDF(mplfile)
    
    if time_resample:
        MPLdat_event.time_resample(timestep,starttime=starttime,endtime=endtime)
    
    if alt_resample:
        MPLdat_event.alt_resample(altrange)
        
    if NRBmask:
        NRBmaskkwargs={'NRBmasktype':NRBmasktype,'NRBthreshold':NRBthresh,'NRBmin':NRBmin,
                       'minalt':NRBminalt,'numprofs':NRBnumprofs,'winsize':NRBwinsize,
                       'inplace':NRBinplace}
        MPLdat_event=mtools.NRB_mask_all(MPLdat_event,**NRBmaskkwargs)


    if doplot:
        plotfilename = '{0}.png'.format(mplfile.split('.')[0])
        plotkwargs={'topplot_limits':NRB_limits,'bottomplot_limits':depol_limits,
                    'hours':hours,'fsize':fsize,'savefilename':plotfilename,'savefilepath':plotsavepath,
                    'SNRmask':SNRmask,'SNRthresh':SNRthresh,'SNRtype':SNRtype,'interpolate':interpolate,
                    'showplot':showplots,'saveplot':saveplots}
        
        doubleplot(MPLdat_event,plotfilename=plotfilename,**plotkwargs)
    
    if docorplot:
        coefplotfilename = '{0}-coefficients.png'.format(mplfile.split('.')[0])
        coefplotkwargs={'toptype':'backscatter','bottomtype':'extinction',
                        'topplot_limits':back_limits,'bottomplot_limits':ext_limits,
                        'hours':hours,'fsize':fsize,'savefilename':coefplotfilename,
                        'savefilepath':plotsavepath,#'SNRmask':SNRmask,'SNRthresh':SNRthresh,
                        'interpolate':interpolate,'SNRtype':SNRtype,'showplot':showplots,'saveplot':saveplots}
            
        doubleplot(MPLdat_event,plotfilename=coefplotfilename,**coefplotkwargs)

    if dolayerplot:
        layerplotfilename = '{0}-layers.png'.format(mplfile.split('.')[0])
        layerplotkwargs={'hours':hours,'fontsize':fsize,'plotfilename':layerplotfilename,
                         #'SNRmask':SNRmask,'SNRthresh':SNRthresh,'SNRtype':SNRtype,
                         'plotfilepath':plotsavepath,'showplot':showplots,'saveplot':saveplots}        
        colormask_plot(MPLdat_event,**layerplotkwargs)
    
if __name__=='__main__': 
    alt_resample=True      
    altrange=np.arange(0.150,5.030,0.030)
    time_resample=False
    timestep='600S'
    starttime=datetime.datetime(2015,07,06,00)
    endtime=datetime.datetime(2015,07,06,18)
    plotsavepath='..\Figures'
    
    NRBmask=False
    NRBthresh=3.0
    NRBmin=0.5
    NRBminalt=0.150
    NRBnumprofs=1
    SNRmask=True
    SNRthresh=0.5
    SNRtype='NRB'
    saveplots=True
    showplots=True
    dolayerplot=True
    docorplot=True
#    hours=['02','04','06','08','10','12','14','168','18','20','22']
#    hours=['03','06','09','12','15','18','21']
    hours=['06','12','18']
    NRB_limits=(0.0,0.5,0.1) 
    depol_limits=(0.0,0.5,0.1)
    back_limits=(0.0,2e-3,5e-4)
    ext_limits=(0.0,1e-1,2e-2)
    interpolate='none'

    os.chdir('K:\Smoke2015\Processed')

    filepath = mtools.get_files('Select Processed MPL file(s)', filetype = ('.h5', '*.h5'))
    
    kwargs={'alt_resample':alt_resample,'altrange':altrange,
            'time_resample':time_resample,'timestep':timestep,'starttime':starttime,'endtime':endtime,
            'SNRmask':SNRmask,'SNRthresh':SNRthresh,'SNRtype':SNRtype,
            'saveplots':saveplots,'showplots':showplots,
            'dolayerplot':dolayerplot,'docorplot':docorplot,
            'NRBmask':NRBmask,'NRBthresh':NRBthresh,'NRBmin':NRBmin,
            'NRBminalt':NRBminalt,'NRBnumprofs':NRBnumprofs,
            'NRB_limits':NRB_limits,'depol_limits':depol_limits,
            'back_limits':back_limits,'ext_limits':ext_limits,
            'interpolate':interpolate,'hours':hours}

    plotall(filepath[0],**kwargs)
    

    
