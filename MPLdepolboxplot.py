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
from matplotlib import cm
import MPLhisttools as mhist
import json

olddir = os.getcwd()

#os.chdir('C:\SigmaMPL\DATA')

os.chdir('K:\Smoke2015\Processed')

filepath = mtools.get_files('Select MPL file', filetype = ('.h5', '*.h5'))

depollist = []
depoldates = []
jsondat={}

for f in filepath:
    MPLfile = mtools.MPL()
    
    MPLfile.fromHDF(f)
    
    altrange = np.arange(0.150,3.03,0.030)
    print 'Resampling {0}'.format(f)
    MPLfile.alt_resample(altrange)
    
    depolraw=MPLfile.depolrat[0]
    tempdate=depolraw.index[0].strftime("%m-%d")
    print 'Filtering data from {0}'.format(tempdate)
    scenefilt=MPLfile.scenepanel[0]['Sub-Type']
    depolfilt=mhist.scenefilter(depolraw,scenefilt,filterterms=['Smoke / Urban','Polluted Dust'])    
    tempfilt=depolfilt.stack().dropna()
    #clean up obvious outliers (non-physical results)
    tempclean=tempfilt[tempfilt>0.0][tempfilt<=1.0]
    jsondat[tempdate]={'median':tempclean.median(),'std':tempclean.std(),'counts':len(tempclean),'freq':(100.0*len(tempclean)/(len(depolraw.stack().dropna())))}
    depoldates.append(tempdate)    
    depollist.append(tempclean)

fignum = 0

font = {'family' : 'serif',
        'weight' : 'medium',
        'size'   : 22}

matplotlib.rc('font', **font)

print 'Generating Figure'

fignum+=1
fig=plt.figure(fignum,figsize=(45,5),dpi=100)
fig.clf()
the_axis=fig.add_subplot(111)
bplot = the_axis.boxplot(depollist, sym = "")
plt.xticks(range(1,(len(depoldates)+1)),depoldates,rotation=45,ha='right')
the_axis.set_ylabel('Volume Depolarization Ratio', fontsize = 24)
#t = the_axis.set_title('Depolarization Boxplot for Aksu Lidar, 150-4000m', fontsize = 24)
the_axis.set_ylim(0.0,0.5)

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

savetime=[depoldates[0],depoldates[-1]]
fig.savefig('../Figures/{0}-{1}-{2}km-depolboxplot.png'.format(savetime[0],savetime[1],altrange[-1]),
            bbox_inches='tight')
fig.canvas.draw()

with open('{0}-{1}-{2}km-depolstats.txt'.format(savetime[0],savetime[1],altrange[-1]),'w') as outfile:
    json.dump(jsondat,outfile,sort_keys=True,indent=4)