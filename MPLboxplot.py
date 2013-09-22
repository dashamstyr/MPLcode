# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 22:47:46 2013

@author: dashamstyr
"""

import MPLtools as mtools
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import slowhist as h2d
from matplotlib import cm


olddir = os.getcwd()

#os.chdir('C:\SigmaMPL\DATA')

os.chdir('C:\SigmaMPL\DATA\Processed')

filepath = mtools.get_files('Select MPL file', filetype = ('.h5', '*.h5'))

depollist = []
depoldates = []

for f in filepath:
    MPLfile = mtools.MPL()
    
    MPLfile.fromHDF(f)
    
    altrange = np.arange(150,4000,30)
    
    MPLfile.alt_resample(altrange)
    
#    MPLfile.time_resample('10M')
    
    copol = MPLfile.NRB[0]
    crosspol = MPLfile.NRB[1]
    
    depoldates.append(copol.index[0].strftime("%m-%d"))
    
    copolvals = np.hstack(copol.values).astype('float32')
    crosspolvals = np.hstack(crosspol.values).astype('float32')
    
    depolMPL = crosspol.values/copol.values
    
    depolvals = depolMPL/(depolMPL+1)
    
    print copol.index[0].strftime("%m-%d")
    mdat = np.ma.masked_array(depolvals.ravel(),np.isnan(depolvals.ravel()))
    print mdat.mean()
    print mdat.std()
    
    depollist.append(depolvals.ravel())

fignum = 0

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 24}

matplotlib.rc('font', **font)

print depollist

fignum+=1
fig=plt.figure(fignum)
fig.clf()
the_axis=fig.add_subplot(111)
bplot = the_axis.boxplot(depollist, sym = "")
plt.xticks(range(1,(len(depoldates)+1)),depoldates)
the_axis.set_ylabel('Volume Depolarization Ratio', fontsize = 28, weight = 'bold')
#t = the_axis.set_title('Depolarization Boxplot for Aksu Lidar, 150-4000m', fontsize = 24)
plt.ylim((0,0.5))

#t.set_y(1.03)

for k in bplot.iterkeys():
    for line in bplot[k]:
        line.set_linewidth(4)

for line in the_axis.xaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(4)

for line in the_axis.yaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(4)
        
[i.set_linewidth(4) for i in the_axis.spines.itervalues()]

#fig.savefig('{0}_{1}-{2}m-copoldepolraw.png'.format(savetime,altrange[0],altrange[-1]))
fig.canvas.draw()