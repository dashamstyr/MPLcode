def fileproc(newdir, altrange = [], timestep = []):
#    clearall()
    import pandas as pan
    import numpy as np
    import os
    import MPLtools as mtools
    from datetime import datetime

    
    #----------------------------------------------------------------------------
    #Uses tools created in MPL_tools to open all files in a folder and resample
    #them to a regular spacing in altitude/date the concatenates them into one
    #pandas dataframe and plots it using LNC_plot
    #July 05, 2012
    #----------------------------------------------------------------------------
    
    os.chdir(newdir)
    
    #newdir = mtools.set_dir('Select Event Folder')
    #
    #os.chdir(newdir)
    #
    #files = os.listdir(newdir)
    
    files = mtools.get_files('Select MPL files', filetype = ('.mpl', '*.mpl'))
    rawfiles = []
    
    #set altitude range and date step sizes
    
#    altrange = np.arange(150,15000,30)#meters
#    timestep = '600S' #seconds
#    starttime = datetime(2013,9,17,0,0,0)
#    endtime = datetime(2013,9,18,0,0,0)
#    timerange1 = [starttime,endtime]
    
    #check to see if each file has been processed before and separate processed
    #files into a new list
    
    for f in files:
        if '.mpl' in f:
            rawfiles.append(f)
    
    #open, altitude resample, and concatenate data and mask files
    
    print rawfiles
    for r in rawfiles:
        MPLdat_temp = mtools.MPL()
        MPLdat_temp.fromMPL(r)
        if len(altrange) > 0:
            MPLdat_temp.alt_resample(altrange)
    
        try:
            MPLdat_event.append(MPLdat_temp)
        except NameError:
            MPLdat_event = MPLdat_temp
       
    #sort by index to make certain data is in order then set date ranges to match
    MPLdat_event.header.sort_index()
    
    for n in range(MPLdat_event.header['numchans'][0]):
        data = MPLdat_event.data[n]
        data = data.sort_index()
    
#    start = data.index[0]
#    end = data.index[-1]
    
    if len(timestep) > 0:
        MPLdat_event.time_resample(timestep)
    
    MPLdat_event.range_cor()    
    MPLdat_event.calculate_NRB(showplots = False)
    MPLdat_event.calculate_depolrat()
    
    olddir = os.getcwd()
    
    try:
        os.chdir('Processed')
    except WindowsError:
        os.makedirs('Processed')
        os.chdir('Processed')
    
    if len(rawfiles) == 1:   
        [path,startfile] = os.path.split(rawfiles[0])
        d_filename = startfile.split('.')[0]+'_proc'
    else:
        [path,startfile] = os.path.split(rawfiles[0])
        [path,endfile] = os.path.split(rawfiles[-1])
        d_filename = startfile.split('.')[0]+'-'+endfile.split('.')[0]+'_proc'
    
    print 'Saving '+d_filename
    MPLdat_event.save_to_HDF(d_filename+'.h5')
    print 'Done'
    os.chdir(olddir)
    

def quickplot(filename, datadir = [], savefigs = False):
    
    import MPLtools as mtools
    import matplotlib.pyplot as plt
    import MPL_plot as mplot
    import os
    import numpy as np
    import pandas as pan
    
#    MPLdat = mtools.MPL()
#    
#    MPLdat.fromHDF(filename[0])
    
    NRBcopol = pan.read_hdf(filename[0],'LNC')
    NRBdepol = pan.read_hdf(filename[0],'MPL')    
    
    if savefigs:    
        if not datadir:
            datadir = mtools.set_dir("Select Folder for Figures to e depositied")
        
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
    
    
def clearall():
    """clear all globals"""
    for uniquevar in [var for var in globals().copy() if var[0] != "_" and var != 'clearall']:
        del globals()[uniquevar]
        
if __name__ == '__main__':
    import MPLtools as mtools
    import os
    import numpy as np
    
    os.chdir('C:\SigmaMPL\DATA')
    procdir = mtools.set_dir('Select folder to process files from')
    altrange = np.arange(150,8000,30)
    timestep = '5Min'
    fileproc(procdir, altrange, timestep)
    
#    plotfile = mtools.get_files('Select File to Plot', filetype = ('.h5', '*.h5'))
#    quickplot(plotfile)
    