def depol_plot(fig, ax, ar, xdata, ydata, data, vrange, fsize = 21):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    
    #set colormap to be the same as 'jet' with the addition of white color for
    #depol ratios set to identically zero because they couldn't be calculated
    cdict = {'red': ((0,1,1),
                     (0.0001, 1, 0),
                     (0.35, 0, 0),
                     (0.66, 1, 1),
                     (0.89,1, 1),
                     (1, 0.5, 0.5)),
         'green': ((0,1,1),
                   (0.0001, 1, 0),
                   (0.125,0, 0),
                   (0.375,1, 1),
                   (0.64,1, 1),
                   (0.91,0,0),
                   (1, 0, 0)),
         'blue': ((0,1,1),
                  (0.0001,1,0.5),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65,0, 0),
                  (1, 0, 0))}
    
    my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,1064)
    
    im = ax.imshow(data, vmin=vrange[0], vmax=vrange[1], cmap = my_cmap)
    forceAspect(ax,ar)
        
    altticks(ax, ydata, fsize = fsize)

    
    ax.set_ylabel('Altitude [m]', fontsize = fsize+4, labelpad = 15)

    for line in ax.yaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(1)
        
    ax.axis('tight')

    return im

def backscatter_plot(fig, ax, ar, xdata, ydata, data, vrange, fsize = 21):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    
    #set colormap to be the same as 'jet' with the addition of white color for
    #depol ratios set to identiacally zero because they couldn't be calculated
    cdict = {'red': ((0,0,0),
                     (0.0099,0,0),
                     (0.01, 1, 0),
                     (0.35, 0, 0),
                     (0.66, 1, 1),
                     (0.89,1, 1),
                     (1, 0.5, 0.5)),
         'green': ((0,0,0),
                   (0.0099,0,0),
                   (0.01, 1, 0),
                   (0.125,0, 0),
                   (0.375,1, 1),
                   (0.64,1, 1),
                   (0.91,0,0),
                   (1, 0, 0)),
         'blue': ((0,0,0),
                  (0.0099,0,0),
                  (0.01,1,0.5),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65,0, 0),
                  (1, 0, 0))}
    
    my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,1064)
    
    im = ax.imshow(data, vmin=vrange[0], vmax=vrange[1], cmap = plt.cm.jet)
    forceAspect(ax,ar)       
    altticks(ax, ydata, fsize = fsize, tcolor = 'w')

    
    ax.set_ylabel('Altitude [m]', fontsize = fsize+4, labelpad = 15)

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
    import matplotlib.pyplot as plt
    from time import strftime
    
    dold = axisdat[0].strftime('%d')
    hold = axisdat[0].strftime('%H')
    tickmarks = []
    ticklabels = []
    n = 0
    l = len(axisdat)

    
    for d in axisdat:
        dtemp = d.strftime('%d')
        if dtemp != dold:
            ticklabels.append(d.strftime('%H\n%b %d'))
            tickmarks.append(n)
        else:
            htemp = d.strftime('%H')
            mtemp = d.strftime('%M')
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
    
    plt.xticks(tickmarks,ticklabels,fontsize = fsize)

    for line in ax.xaxis.get_ticklines():
        line.set_color(tcolor)
        line.set_markersize(10)
        line.set_markeredgewidth(2)

def altticks(ax, axisdat, numticks = 5, fsize = 21, tcolor = 'k'):
    import matplotlib.pyplot as plt

    numpoints = len(axisdat)
    step = numpoints//numticks
    tickmarks = range(0,numpoints,step)
    ticklabels = [str(int(t)) for t in axisdat[::step]]

    plt.yticks(tickmarks,ticklabels, fontsize = fsize)

    for line in ax.yaxis.get_ticklines():
        line.set_color(tcolor)
        line.set_markersize(10)
        line.set_markeredgewidth(3)

def vertprof(df, altrange, exact_times, plot_type = 'line', zeromask = False, 
             savefig = False, filename = 'tempfig.png'):    
    import os, sys
    import numpy as np
    import datetime as dt
    import bisect
    import matplotlib.pyplot as plt
    import matplotlib.colors as clr
    import pandas as pan
    import LNC_tools as LNC
    import LNC_plot as LNCplt
    
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
        
    
if __name__ == '__main__':
    import pandas as pan
    import os
    import MPLtools as MPL
    import matplotlib.pyplot as plt
    import numpy as np

    os.chdir('C:\Program Files (x86)\SigmaMPL\DATA')

    filepath = mtools.get_files('Select MPL file', filetype = ('.h5', '*.h5'))
    
    MPLfile = mtools.MPL()
    
    MPLfile.fromHDF(filepath[0])
    
    altrange = np.arange(150,20000,10)
    
    MPLfile.alt_resample(altrange)
    
    #create figure and plot image of depolarization ratios
    #create figure and plot image of depolarization ratios
    fsize = 18 #baseline font size
    ar = 2.0  #aspect ratio
    figheight = 12 #inches
    
    plt.rc('font', family='serif', size=fsize)
    
    
    fig = plt.figure()
    
    h_set = range(1,25)
    h_set = map(str,h_set)
    
    copol = MPLfile.data[0]
    crosspol = MPLfile.data[1]
    
    depolMPL = crosspol.values/copol.values
    
    depolvals = depolMPL/(depolMPL+1)
    
    depol = pan.DataFrame(depolvals,index = copol.index, columns = copol.columns)
    #depol = imtools.blur_image(depol,(7,7), kernel ='Flat')
    
    datetime = copol.index
    alt = copol.columns
    
    print 'Generating Figure'
    
    ax1 = fig.add_subplot(2,1,1)
    im1 = backscatter_plot(fig, ax1, ar, datetime,alt[::-1],copol.T[::-1], fsize = fsize)
    cbar1 = fig.colorbar(im1, orientation = 'vertical', aspect = 6)
    cbar1.ax.tick_params(labelsize = fsize)
    dateticks(ax1, datetime, hours = h_set, fsize = fsize, tcolor = 'w')
    ax1.set_xticklabels([])
    t1 = ax1.set_title('Attenuated Backscatter', fontsize = fsize+10)
    t1.set_y(1.03)
            
    ax2 = fig.add_subplot(2,1,2)
    im2 = depol_plot(fig, ax2, ar, datetime,alt[::-1],depol.T[::-1], fsize = fsize)
    cbar2 = fig.colorbar(im2, orientation = 'vertical', aspect = 6)
    cbar2.set_ticks(np.arange(0,1.1,0.1))
    cbar2.set_ticklabels(np.arange(0,1.1,0.1))
    cbar2.ax.tick_params(labelsize = fsize)
    #set axis ranges and tickmarks based on data ranges
    dateticks(ax2, datetime, hours = h_set, fsize = fsize)
    ax2.set_xlabel('Time [Local]',fontsize = fsize+4)
    fig.autofmt_xdate()
    t2 = ax2.set_title('Linear Depolarization Ratio',fontsize = fsize+10)
    t2.set_y(1.03)    
    
    ##plt.savefig(savetitle,dpi = 100, edgecolor = 'b', bbox_inches = 'tight')
    fig.set_size_inches(figheight*ar,figheight) 
    plt.savefig(d_filename+'.png')
    print 'Done'
    
    plt.show()
    os.chdir(olddir)