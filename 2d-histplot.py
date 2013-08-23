# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:10:30 2013

@author: dashamstyr
"""

import MPLtools as mtools
import os
import numpy as np
import matplotlib.pyplot as plt
import slowhist as h2d
from matplotlib import cm


olddir = os.getcwd()

#os.chdir('C:\SigmaMPL\DATA')

os.chdir('C:\SigmaMPL\DATA\Processed')

filepath = mtools.get_files('Select MPL file', filetype = ('.h5', '*.h5'))

MPLfile = mtools.MPL()

MPLfile.fromHDF(filepath[0])

altrange = np.arange(150,8000,30)

MPLfile.alt_resample(altrange)

copol = MPLfile.NRB[0]
crosspol = MPLfile.NRB[1]

copolvals = np.hstack(copol.values).astype('float32')
crosspolvals = np.hstack(crosspol.values).astype('float32')

depolMPL = crosspol.values/copol.values

depolvals = depolMPL/(depolMPL+1)

numbins = 100
depolmin = 0.0
depolmax = 0.5
copolmin = 0.0
copolmax = 3e-2

copolhist=h2d.fullhist(copolvals,numbins,copolmin,copolmax,-9999.,-8888.)
depolhist=h2d.fullhist(np.hstack(depolvals),numbins,depolmin,depolmax,-9999.,-8888.)

altOut = h2d.althist(depolvals,altrange,numbins,(depolmin,depolmax))

copolOut=h2d.hist2D(copolhist['fullbins'],depolhist['fullbins'],copolhist['numBins'],depolhist['numBins'])

altcounts=altOut['coverage']
copolcounts = copolOut['coverage']

altcounts[altcounts < 1] = 1
copolcounts[copolcounts < 1] = 1

altlogcounts=np.log10(altcounts)
copollogcounts=np.log10(copolcounts)
                  
try:
    os.chdir('../HistPlots')
except WindowsError:
    os.makedirs('../HistPlots')
    os.chdir('../HistPlots')

startdate = copol.index[0].strftime("%Y-%m-%d")
enddate = copol.index[-1].strftime("%Y-%m-%d")

starttime = copol.index[0].strftime("%H")
endtime = copol.index[-1].strftime("%H")

if startdate == enddate:
    if starttime == endtime:
        savetime = startdate+'_'+starttime
    else:
        savetime = startdate+'_'+starttime+'-'+endtime
else:
    savetime = startdate+'-'+enddate

fignum = 0

#fignum+=1
#fig=plt.figure(fignum)
#fig.clf()
#the_axis=fig.add_subplot(111)
#the_axis.plot(depolvals,copolvals,'b+')
#the_axis.set_xlabel('depolvals')
#the_axis.set_ylabel('copolvals')
#the_axis.set_title('raw scatterplot')
#fig.savefig('{0}_{1}-{2}m-copoldepolraw.png'.format(savetime,altrange[0],altrange[-1]))
#fig.canvas.draw()

cmap=cm.bone
cmap.set_over('r')
cmap.set_under('b')

fignum+=1
fig=plt.figure(fignum)
fig.clf()
axis=fig.add_subplot(111)
im=axis.pcolormesh(depolhist['centers'],altrange,altlogcounts.T, cmap = cmap)
cb=plt.colorbar(im,extend='both')
title="2-D Histogram: Altitude vs. Depol. Ratio"
colorbar="log10(counts)"
the_label=cb.ax.set_ylabel(colorbar,rotation=270)
axis.set_xlabel('Volume Depolarization Ratio')
axis.set_ylabel('Altitude [m]')
axis.set_title(title)
fig.savefig('{0}_{1}-{2}m-althist.png'.format(savetime,altrange[0],altrange[-1]))
fig.canvas.draw()

fignum+=1
fig=plt.figure(fignum)
fig.clf()
axis=fig.add_subplot(111)
im=axis.pcolormesh(depolhist['centers'],copolhist['centers'],copollogcounts, cmap = cmap)
cb=plt.colorbar(im,extend='both')
title="2-D Histogram: Backscatter vs. Depol. Ratio"
colorbar="log10(counts)"
the_label=cb.ax.set_ylabel(colorbar,rotation=270)
axis.set_xlabel('Volume Depolarization Ratio')
axis.set_ylabel('Attenuated Backscatter $[km^{-1} sr^{-1}]$')
axis.set_title(title)
fig.savefig('{0}_{1}-{2}m-copoldepol.png'.format(savetime,altrange[0],altrange[-1]))
fig.canvas.draw()

#fignum+=1
#fig = plt.figure(fignum)
#fig.clf()
#axis = fig.add_subplot(111)
#hist = axis.hist(depolvals, bins = 100)
#axis.set_xlabel('Volume Depolaization Ratio')
#fig.savefig('{0}_{1}-{2}m-1Ddepolhist.png'.format(savetime,altrange[0],altrange[-1]))
#fig.canvas.draw()
#plt.show()

os.chdir(olddir)