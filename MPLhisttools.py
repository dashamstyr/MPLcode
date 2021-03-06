
import os,sys,site
#home=os.environ['homepath']
#site.addsitedir('{0}\\Dropbox\\Python_Scripts\\GIT_Repos\\'.format(home))

import numpy as np
from scipy.stats import cumfreq
import time, datetime
from math import isnan
import pandas as pan
from matplotlib import pyplot as plt
from matplotlib import cm, ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from copy import deepcopy
from scipy import signal

import MPLtools as mtools
import MPLprocesstools as mproc
import MPLfileproc as mfile
import MPLplot as mplot
import LNC_tools as ltools

def progress(datanum,tot_loops):
    """
       print a message about percent complete
       input:      datanum: current index (int)
       tot_loops:  index when loop finishes
    """   
    the_frac=np.int(np.float(datanum)/tot_loops*100.0)
    sys.stdout.write("\rpercent complete: {0}%".format(the_frac))
    sys.stdout.write("{0} files out of {1}".format(datanum,tot_loops))
    sys.stdout.flush()

def filegetter(filedir=None,savedir=None,filetype='.mpl',**kwargs):
    """
    takes a top-level directory and finds all relevant files 
    returns a list of filenames to iterate through
    
    inputs:
    filedir = top-level directory.  If emplty, interactive window is used
    filetype = file extension to look for (.mpl or .h5), defaults to .mpl
    """    
    altitudes = kwargs.get('altitudes',None)
    starttime = kwargs.get('starttime',datetime.datetime(1970,1,1,0))
    endtime = kwargs.get('endtime',datetime.datetime.now())
    timestep = kwargs.get('timestep',None)
    verbose = kwargs.get('verbose',False)
    SNRmask=kwargs.get('SNRmask',True)
    SNRthreshold=kwargs.get('SNRthreshold',1.0)
    bg_alt=kwargs.get('bg_alt',None)
    datatype=kwargs.get('datatype','NRB')
    NRBmask=kwargs.get('NRBmask',True)
    NRBthreshold=kwargs.get('NRBthreshold',5.0)
    NRBmin=kwargs.get('NRBmin',0.5)
    minalt=kwargs.get('minalt',0.150)
    dolayers=kwargs.get('dolayers',True)
    docorrection=kwargs.get('docorrection',False)
    winsize=kwargs.get('winsize',5)
    topickle=kwargs.get('topickle',True)
    picklefilename=kwargs.get('picklefilename','test.p')
    
    #findalllayers kwargs
    molthresh=kwargs.get('molthresh',1.0)
    layernoisethresh=kwargs.get('layernoisethresh',1.0)      
    bg_alt=kwargs.get('bg_alt',None)
    datatype=kwargs.get('datatype','NRB')
    winsize=kwargs.get('winsize',5)
    wavelet=kwargs.get('wavelet',signal.ricker)
    noisethresh=kwargs.get('noisethresh',0.4)
    cloudthresh=kwargs.get('cloudthresh',(1.0,0.20))
    CWTwidth=kwargs.get('CWTwidth',2)
    minwidth=kwargs.get('minwidth',4)
    layerCWTrange=kwargs.get('layerCWTrange',np.arange(2,5))
    PBLwavelet=kwargs.get('PBLwavelet',mproc.dog)
    PBLCWTrange=kwargs.get('PBLCWTrange',np.arange(5,15))
    PBLwidth=kwargs.get('PBLwidth',7)
    sigma0=kwargs.get('sigma0',0.1)
    depolsigma0=kwargs.get('depolsigma0',0.05)
    waterthresh=kwargs.get('waterthresh',0.10)
    icethresh=kwargs.get('icethresh',0.35)
    smokethresh=kwargs.get('smokethresh',0.05)
    dustthresh=kwargs.get('dustthresh',0.20)
    maxaeroalt=kwargs.get('maxaeroalt',10.0)
    
    #correction kwargs
    refalt=kwargs.get('refalt',None)
    calrange=kwargs.get('calrange',None)
    method=kwargs.get('method','klett2')
    lrat=kwargs.get('lrat',None)
    mode=kwargs.get('mode','copol')
    
    doprogress=kwargs.get('doproc',True)
    recalc=kwargs.get('recalc',False)
    
    if filedir is None:
        filedir = mtools.set_dir('Sepect top-level folder to search')
    
    if savedir is None:
        savedir = mtools.set_dir('Sepect top-level folder to save to')
        
    outlist=[]
    picklelist=[]
    n=0
    for root, folder, filelist in os.walk(filedir):
        for filename in filelist:
            if filename.endswith(filetype):
                if doprogress:
                    progress(n,len(filelist))
                    n+=1
                tempMPL=mtools.MPL()
                if filetype=='.mpl':   
                    if starttime<=mfile.MPLtodatetime(filename)<=endtime:
                        fileinfo=os.stat(os.path.join(root,filename))
                        if fileinfo.st_size==0:
                            continue
                        else:
                            savefilename='{0}_proc.h5'.format(filename.split('.')[0])
                            if recalc:
                                tempMPL.fromMPL(os.path.join(root,filename))
                                procflag=False
                            elif savefilename in os.listdir(savedir):
                                procflag=True
                                tempMPL.fromHDF(os.path.join(savedir,savefilename))
                            else:
                                tempMPL.fromMPL(os.path.join(root,filename))
                                procflag=False
                        if verbose:
                            print filename
                    else:
                        continue
                elif filetype=='.h5':
                    procflag=True
                    if starttime<=mfile.H5todatetime(filename)<=endtime:
                        fileinfo=os.stat(os.path.join(root,filename))
                        if fileinfo.st_size==0:
                            continue
                        else:
                            tempMPL.fromHDF(os.path.join(root,filename))
                        if verbose:
                            print filename
                    else:
                        continue 
                if procflag and not recalc:
                    if topickle:
                        picklelist.append(tempMPL)
                    outlist.append(os.path.join(savedir,savefilename))
                        
                else:
                    if timestep:
                        tempMPL = tempMPL.time_resample(timestep=timestep,verbose=verbose)
                    if altitudes.any():
                        tempMPL = tempMPL.alt_resample(altitudes,verbose=verbose)
                
                    tempMPL.calc_all()
                    if dolayers:
                        layerkwargs={'timestep':timestep,'bg_alt':bg_alt,'datatype':datatype,
                                     'molthresh':molthresh,'winsize':winsize,'layernoisethresh':layernoisethresh,
                                     'wavelet':wavelet,'noisethresh':noisethresh,'cloudthresh':cloudthresh,
                                     'CWTwidth':CWTwidth,'minwidth':minwidth,'layerCWTrange':layerCWTrange,
                                     'PBLwavelet':PBLwavelet,'PBLCWTrange':PBLCWTrange,'PBLwidth':PBLwidth,'sigma0':sigma0,
                                     'waterthresh':waterthresh,'icethresh':icethresh,'smokethresh':smokethresh,
                                     'dustthresh':dustthresh,'sigma0':sigma0,'depolsigma0':depolsigma0,'maxaeroalt':maxaeroalt}
                        layerdict=mproc.findalllayers(mplin=tempMPL,**layerkwargs)
                        tempMPL=mproc.scenemaker(layerdict)
                    
                    if docorrection:
                        corkwargs={'refalt':refalt,'calrange':calrange,'method':method,'lrat':lrat,'mode':mode}
                        
                        tempMPL=mproc.basiccorrection(tempMPL,**corkwargs)
                        
                    if NRBmask:
                        tempMPL=mtools.NRB_mask_all(tempMPL,NRBthreshold=NRBthreshold,
                                                   NRBmin=NRBmin,minalt=minalt,
                                                   winsize=winsize)
                    if SNRmask:
                        tempMPL=mtools.SNR_mask_all(tempMPL,SNRthreshold=SNRthreshold,
                                                   bg_alt=bg_alt,datatype=datatype)                            
                        
                    tempMPL.save_to_HDF(os.path.join(savedir,savefilename))
                    if topickle:
                        picklelist.append(tempMPL)
                    outlist.append(os.path.join(savedir,savefilename))
            
    if topickle:
        pickle.dump(picklelist,open(picklefilename,'wb'))
    
    return outlist
        
def dataextractor(datatype='NRB',**kwargs):
    """
    takes a file name, or a list of MPL files
    opens the file and extracts the MPL data into a dataframe
    and filters data based on kwargs returns an MPL class object or dataframe
    
    inputs:
    datatype = type of data to extract (defaults to NRB)
    
    kwargs:
    loadfiletype
    
    """
    loadfiletype=kwargs.get('loadfiletype','list')
    loadfilename=kwargs.get('loadfilename',None)
    MPLlist=kwargs.get('MPLlist',None)
    toHDF=kwargs.get('toHDF',False)
    HDFfile=kwargs.get('HDFfile','test.h5')
    doprogress=kwargs.get('doprogress',True)
    dflist=[]

    if loadfiletype=='HDF5':
        try:
            mpltemp=mtools.MPL()
            mpltemp.fromHDF(loadfilename)
            if datatype=='Raw':
                dfout=mpltemp.data[0]
            elif datatype=='NRB':
                dfout=mpltemp.NRB[0]
            elif datatype=='Depolrat':
                dfout=mpltemp.depolrat[0]
            elif datatype=='SNR':
                dfout=mpltemp.SNR['NRB'][0]
            elif datatype=='DepolSNR':
                dfout=mpltemp.SNR['depolrat'][0]
            elif datatype=='Type':
                dfout=mpltemp.scenepanel[0]['Type']
            elif datatype=='Sub-Type':
                dfout=mpltemp.scenepanel[0]['Sub-Type']
            elif datatype=='Delta':
                dfout=mpltemp.scenepanel[0]['Delta']
            
        except IOError:
            print "Could not find file named {0}".format(loadfilename)
    elif loadfiletype=='LNCHDF':
        try:
            LNCtemp=ltools.LNC()
            LNCtemp.fromHDF(loadfilename)
            if datatype=='BR':
                dfout=LNCtemp.BR
            elif datatype=='PR':
                dfout=LNCtemp.PR
            elif datatype=='SNR':
                dfout=LNCtemp.SNR['BR']
            elif datatype=='Type':
                dfout=LNCtemp.scenepanel['Type']
            elif datatype=='Sub-Type':
                dfout=LNCtemp.scenepanel['Sub-Type']
            elif datatype=='Delta':
                dfout=LNCtemp.scenepanel['Delta']
        except IOError:
            print "Could not find file named {0}".format(loadfilename)
    else:            
        if loadfiletype=='pickle':
            try:
                with open(loadfilename,'rb') as pf:
                    MPLlist=pickle.load(pf)
            except IOError:
                print "Could not find file named {0}".format(loadfilename)
        if MPLlist is not None:
            for MPLfile in MPLlist:        
                if datatype=='Raw':
                    tempdf=MPLfile.data[0]
                elif datatype=='NRB':
                    tempdf=MPLfile.NRB[0]
                elif datatype=='Depolrat':
                    tempdf=MPLfile.depolrat[0]
                elif datatype=='SNR':
                    tempdf=MPLfile.SNR['NRB'][0]
                elif datatype=='DepolSNR':
                    tempdf=MPLfile.SNR['depolrat'][0]
                elif datatype=='Type':
                    tempdf=MPLfile.scenepanel[0]['Type']
                elif datatype=='Sub-Type':
                    tempdf=MPLfile.scenepanel[0]['Sub-Type']
                elif datatype=='Delta':
                    dfout=mpltemp.scenepanel[0]['Delta']
                dflist.append(tempdf)
            
            dfout=pan.concat(dflist)
        else:
            print 'Input file type {0} not recognized!'.format(loadfiletype)
            return dflist
            
    if toHDF:
        store=pan.HDFStore(HDFfile)
        store[datatype]=dfout
        store.close()
    
    return dfout

def scenefilter(dfin,dffilt,**kwargs):
    filtmode=kwargs.get('filtmode','Type')
    filterterms=kwargs.get('filterterms',None)
        
    if filtmode=='Type':    
        dftemp=dffilt.isin(filterterms)
        dftemp.replace(False,np.nan,inplace=True)
        dfout=dfin*dftemp
    elif filtmode=='LT':
        dfout=dfin[dffilt<=filterterms[0]]
        dfout.replace(False,np.nan,inplace=True)
    elif filtmode=='GT':
        dfout=dfin[dffilt>=filterterms[0]]
        dfout.replace(False,np.nan,inplace=True)
    elif filtmode=='Range':
        dfout=dfin[dffilt>=filterterms[0]]
        dfout=dfout[dfout<=filterterms[1]]
        dfout.replace(False,np.nan,inplace=True)
       
    return dfout

#def SNRfilter(dfin,dfSNR,thresh=1.0):
#    if inplace:
#        dftemp=dfin
#    else:
#        dftemp=deepcopy(dfin)
#    
#    dfout=dfin[dfSNR>thresh]
#    dfout.replace(False,np.nan,inplace=True)
#        
#    return dfout

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
    
def histplot1D(datain,**kwargs):
    
    datatype=kwargs.get('datatype','df')
    
    if datatype is 'df':
        histvals=datain.values
        binrange=kwargs.get('binrange',[datain.min(),datain.max()])
    elif datatype is 'histdict':
        counts=datain['counts']
        edges=datain['edges']
        centers=datain['centers']
        step=np.diff(edges)
        
    numbins=kwargs.get('numbins',100)  
    missinglowval=kwargs.get('missinghighval',-99999)
    missinghighval=kwargs.get('missinglowval',99999)
    normalize=kwargs.get('normalize',True)
    cumulative=kwargs.get('cumulative',False)
    
    doplot=kwargs.get('doplot',True)  
    showplot=kwargs.get('showplot',False)
    saveplot=kwargs.get('saveplot',False)
    plotfilename=kwargs.get('plotfilename','1Dhist_test.png')
    fsize=kwargs.get('fsize',32) #baseline font size
    ar=kwargs.get('ar',1.0)  #aspect ratio
    figheight=kwargs.get('figheight',12) #inches 
    dpi=kwargs.get('dpi',100)      
    fignum=kwargs.get('fignum',0)
    xlog=kwargs.get('xlog',False)
    ylog=kwargs.get('ylog',False)
    xlimits=kwargs.get('xlimits',None)
    ylimits=kwargs.get('ylimits',None)
    xlabel=kwargs.get('xlabel',None)
    if ylog:
        ylabel=kwargs.get('ylabel','Log (Counts)')
    else:
        ylabel=kwargs.get('ylabel','Counts')
    
    if datatype is 'histdict':
        dictout=datain
    else:
        if not cumulative:
            counts,edges=np.histogram(datain.values,numbins,range=binrange,normed=normalize) 
            step=np.diff(edges)
            centers=edges[:-1]+step*0.5            
            dictout={'counts':counts,'centers':centers,'edges':edges}
        else:            
            counts,lowlim,barwidths,extrapoints=cumfreq(datain.values,numbins=numbins,defaultreallimits=binrange)
            if normalize:
                totcounts=((datain.values>missinglowval)&(datain.values<missinghighval)).sum()
                counts=counts/totcounts
            step=(binrange[1]-binrange[0])/numbins
            binvals=np.arange(binrange[0],binrange[1],step)
            centers=[v+step*0.5 for v in binvals]
            edges=np.hstack((binvals,binvals[-1]+step))
            dictout={'counts':counts,'centers':centers,'edges':edges}
    if doplot:
        plt.rc('font', family='serif', size=fsize)
        fig1=plt.figure(fignum)
        if ar:
            fig1.set_size_inches(figheight*ar,figheight)
        ax1=fig1.add_subplot(111)
        barwidths=step
        if xlog:
            logplot=True
            plt.xscale('log')
        else:
            logplot=False
        
        if ylog:
            logplot=True
            plt.yscale('log')
        else:
            logplot=False
            ax1.set_aspect(ar)
            
        ax1.bar(centers,counts,width=barwidths,align='center',log=logplot)  
        

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if ylimits:
            plt.ylim(ylimits)
        if xlimits:
            plt.xlim(xlimits)
        if saveplot:
            fig1.canvas.print_figure(plotfilename,dpi = dpi, edgecolor = 'b', bbox_inches = 'tight') 
        if showplot:
            fig1.canvas.draw()
    return dictout

def histplot2D(dfx,dfy,**kwargs):
    
    xbinrange=kwargs.get('xbinrange',None)
    ybinrange=kwargs.get('ybinrange',None)
    numbins=kwargs.get('numbins',[100,100])  
    missinglowval=kwargs.get('missinghighval',-99999)
    missinghighval=kwargs.get('missinglowval',99999)
    normalize=kwargs.get('normalize',True)
    
    doplot=kwargs.get('doplot',True)    
    saveplot=kwargs.get('savefile',False)
    plotfilename=kwargs.get('plotfilename','2Dhist_test.png')
    fsize=kwargs.get('fsize',32) #baseline font size
    ar=kwargs.get('ar',1.0)  #aspect ratio
    figheight=kwargs.get('figheight',12) #inches
    xlabel=kwargs.get('xlabel',None)
    ylabel=kwargs.get('ylabel',None)        
    fignum=kwargs.get('fignum',0)
    cmap=kwargs.get('cmap',cm.bone)
    undercolor=kwargs.get('undercolor','r')
    overcolor=kwargs.get('overcolor','b')
    logcolor=kwargs.get('logcolor',True)
    
    if len(np.shape(dfx))==1:
        xvalues=dfxdropna().values
    else:
        xvalues=dfx.stack().dropna().values
        
    if len(np.shape(dfy))==1:
        yvalues=dfydropna().values
    else:
        yvalues=dfy.stack().dropna().values
    
    if xbinrange is None:
        xbinrange=[min(xvalues),max(xvalues)]
    
    if ybinrange is None:
        ybinrange=[min(yvalues),max(yvalues)]   
        
    xycounts,xedges,yedges=np.histogram2d(xvalues,yvalues,numbins,range=[xbinrange,ybinrange])
    
    xstep=np.diff(xedges)
    ystep=np.diff(yedges)
    xcenters=xedges[:-1]+xstep*0.5
    ycenters=yedges[:-1]+ystep*0.5
    
    if doplot:

        if logcolor:
            plotcounts=np.log10(xycounts)
            colorbarlabel='$log_(10) Counts$'
        else:
            colorbarlabel='Counts'
        xyextent=[xedges[0],xedges[-1],yedges[0],yedges[-1]]
        if ar:
            ar=ar*(xedges[-1]-xedges[0])/(yedges[-1]-yedges[0])
        else:
            ar=(xedges[-1]-xedges[0])/(yedges[-1]-yedges[0])
        plt.rc('font', family='serif', size=fsize)
        fig=plt.figure(fignum)
        fig.clf()
        axis=fig.add_subplot(111)
        cmap.set_over(undercolor)
        cmap.set_under(overcolor)
        im=axis.imshow(plotcounts.T,origin='low',extent=xyextent,interpolation='none',
                       aspect=ar, cmap = cmap)
        axis.xaxis.set_major_locator(ticker.MaxNLocator(nbins = 5, prune='lower'))
        cb=plt.colorbar(im,extend='both')
        cb.ax.set_ylabel(colorbarlabel,rotation=270,labelpad=50)
        
        for line in axis.yaxis.get_ticklines():
            line.set_color('w')
            line.set_markersize(10)
            line.set_markeredgewidth(3)    
    
        for line in axis.xaxis.get_ticklines():
            line.set_color('w')
            line.set_markersize(10)
            line.set_markeredgewidth(3)
        
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
    #    t2 = axis.set_title(title)
        fig.set_size_inches(figheight*ar,figheight) 
#        fig.tight_layout()
    #    t2.set_y(1.02)
        fig.canvas.draw()
        if saveplot:
            fig.savefig(plotfilename)
    
    dictout={'counts':xycounts,'xedges':xedges,'yedges':yedges,'xcenters':xcenters,
             'ycenters':ycenters}
    return dictout

def althistplot(dfin,**kwargs):
        
    databinrange=kwargs.get('databinrange',[min(dfin.min()),max(dfin.max())])
    altbinrange=kwargs.get('altbinrange',[dfin.columns.min(),dfin.columns.max()])
    numdatbins=kwargs.get('numdatbins',100)  
    numaltbins=kwargs.get('numaltbins',len(dfin.columns))
    
    doplot=kwargs.get('doplot',True)    
    saveplot=kwargs.get('savefile',False)
    plotfilename=kwargs.get('plotfilename','2dalthist_test.png')
    fsize=kwargs.get('fsize',32) #baseline font size
    ar=kwargs.get('ar',None)  #aspect ratio
    figheight=kwargs.get('figheight',12) #inches
    xlabel=kwargs.get('xlabel',None)
    ylabel=kwargs.get('ylabel','Altitude [m]')        
    fignum=kwargs.get('fignum',0)
    cmap=kwargs.get('cmap',cm.bone)
    undercolor=kwargs.get('undercolor','r')
    overcolor=kwargs.get('overcolor','b')
    logcolor=kwargs.get('logcolor',True)
    
    altgrid,dtgrid=np.meshgrid(dfin.columns.values,range(len(dfin.index)))
    datavals=np.ravel(dfin.values)
    altvals=np.ravel(altgrid)
    altcounts,datedges,altedges=np.histogram2d(datavals,altvals,bins=[numdatbins,numaltbins],
                                             range=[databinrange,altbinrange])
    altstep=np.diff(altedges)
    datstep=np.diff(datedges)
    altcenters=altedges[:-1]+altstep*0.5
    datcenters=datedges[:-1]+datstep*0.5
    
    if doplot:
        if logcolor:
            plotcounts=np.log10(altcounts)
            colorbarlabel='$log_(10) Counts$'
        else:
            plotcounts=altcounts
            colorbarlabel='Counts'        
        if ar:
            ar=ar*(datedges[-1]-datedges[0])/(altedges[-1]-altedges[0])
        else:
            ar=(datedges[-1]-datedges[0])/(altedges[-1]-altedges[0])
        plotextent=[datedges[0],datedges[-1],altedges[0],altedges[-1]]
        plt.rc('font', family='serif', size=fsize)
        fignum+=1
        fig=plt.figure(fignum)
        fig.clf()
        axis=fig.add_subplot(111)
        cmap.set_over(undercolor)
        cmap.set_under(overcolor)
        im=axis.imshow(plotcounts.T,origin='low',extent=plotextent,interpolation='none',
                       aspect=ar, cmap = cmap)
        axis.xaxis.set_major_locator(ticker.MaxNLocator(nbins = 5, prune='lower'))
        cb=plt.colorbar(im,extend='both')
        cb.ax.set_ylabel(colorbarlabel,rotation=270,labelpad=50)
        
        for line in axis.yaxis.get_ticklines():
            line.set_color('w')
            line.set_markersize(10)
            line.set_markeredgewidth(3)    
    
        for line in axis.xaxis.get_ticklines():
            line.set_color('w')
            line.set_markersize(10)
            line.set_markeredgewidth(3)
        if xlabel:
            axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        fig.tight_layout()
        
        
        fig.canvas.draw()
        if saveplot:
            fig.savefig(plotfilename)
    
    dictout={'counts':altcounts,'datedges':datedges,'altedges':altedges,'datcenters':datcenters,
             'altcenters':altcenters}
    return dictout

    
def altticks(ax, axisdat, numticks = 5, fsize = 21, tcolor = 'k'):

    numpoints = len(axisdat)
    step = numpoints//numticks
    tickmarks = axisdat[::step]
    ticklabels = [str(int(t)) for t in axisdat[::step]]

    plt.yticks(tickmarks,ticklabels, fontsize = fsize)

    for line in ax.yaxis.get_ticklines():
        line.set_color(tcolor)
        line.set_markersize(10)
        line.set_markeredgewidth(3)
        




if __name__=="__main__":
    

    topdir='C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\DATA'
    savefile='histpickle_local.p'
    os.chdir(topdir)
#    picklefile=mtools.get_files('Select pickle file to histogram',filetype=('.p','*.p'))[0]
    startdate=datetime.datetime(2013,1,1,00)
    enddate=datetime.datetime(2014,11,4,1)
    altitudes=np.arange(0.150,15.000,0.030)
    timestep='600S'
    mplfiles=filegetter(filedir=topdir,altitudes=altitudes,starttime=startdate,
                        endtime=enddate,timestep=timestep,topickle=True,
                        picklefile=savefile,verbose=True)
    
#    with open(savefile,'rb') as sfile:
#        allfiles=pickle.load(sfile)
    
    
    datatypes={'NRB':['NRBhistdat_sigma.h5',(0.0,5.0),'NRBhistplot_sigma.png'],
                      'Depolrat':['Depolhistdat_sigma.h5',(0,1.0),'Depolhistplot_sigma.png'],
                      'SNR':['SNRhistdat_sigma.h5',(0,10.0),'SNRhistplot_sigma.png']}
    alldict={}
    n=0
    for dtype,dspecs in datatypes.iteritems():
        df=dataextractor(datatype=dtype,loadfiletype='pickle',loadfilename=picklefile,
                         toHDF=True,HDFfile=os.path.join(topdir,dspecs[0]))
        dfplot=df[df>0]
        dfplot.fillna(value=-99999,inplace=True)    
        
        alldict[dtype]=histplot1D(df,numbins=100,binrange=dspecs[1],xlabel=dtype,
                                    fignum=n,cumulative=False,normalize=True,
                                    saveplot=True,plotfilename=os.path.join(topdir,dspecs[2]))
        n+=1


    os.chdir(topdir)
    mplfile=mtools.get_files('blah',filetype = ('.h5','*.h5'))
    MPLtest=mtools.MPL()
    MPLtest.fromHDF(mplfile[0])
    MPLtest.calc_all()
    NRB=MPLtest.NRB[0]
    NRB.fillna(-99999,inplace=True)
    depol=MPLtest.depolrat[0]
    depol.fillna(-99999,inplace=True)
    hist1,NRBedges,altedges=althistplot(NRB,numdatbins=500,numaltbins=len(NRB.columns),databinrange=[0.0,10.0],doplot=True)
