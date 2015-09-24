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

NRBlist = []
NRBdates = []
jsondat={}

for f in filepath:
    MPLfile = mtools.MPL()
    
    MPLfile.fromHDF(f)
    
    altrange = np.arange(0.150,3.03,0.030)
    print 'Resampling {0}'.format(f)
    MPLfile.alt_resample(altrange)
    
    NRBraw=MPLfile.NRB[0]
    tempdate=NRBraw.index[0].strftime("%m-%d")
    print 'Filtering data from {0}'.format(tempdate)
    scenefilt=MPLfile.scenepanel[0]['Sub-Type']
    NRBfilt=mhist.scenefilter(NRBraw,scenefilt,filterterms=['Smoke / Urban','Polluted Dust'])    
    tempfilt=NRBfilt.stack().dropna()
    #clean up obvious outliers (non-physical results)
    tempclean=tempfilt[tempfilt>0.0]
    jsondat[tempdate]={'median':tempclean.median(),'std':tempclean.std(),'counts':len(tempclean),'freq':(100.0*len(tempclean)/(len(NRBraw.stack().dropna())))}
    NRBdates.append(tempdate)    
    NRBlist.append(tempclean)

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
bplot = the_axis.boxplot(NRBlist, sym = "")
plt.xticks(range(1,(len(NRBdates)+1)),NRBdates,rotation=45,ha='right')
the_axis.set_ylabel('NRB [${counts*km^{2}}/{\mu s*\mu J}$]', fontsize = 24)
#t = the_axis.set_title('Depolarization Boxplot for Aksu Lidar, 150-4000m', fontsize = 24)
the_axis.set_ylim(0.0,1.0)

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

savetime=[NRBdates[0],NRBdates[-1]]
fig.savefig('../Figures/{0}-{1}-{2}km-NRBboxplot.png'.format(savetime[0],savetime[1],altrange[-1]),
            bbox_inches='tight')
fig.canvas.draw()

with open('{0}-{1}-{2}km-NRBstats.txt'.format(savetime[0],savetime[1],altrange[-1]),'w') as outfile:
    json.dump(jsondat,outfile,sort_keys=True,indent=4)