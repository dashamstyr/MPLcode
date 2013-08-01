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

os.chdir('C:\SigmaMPL\DATA')

filepath = mtools.get_files('Select MPL file', filetype = ('.h5', '*.h5'))

MPLfile = mtools.MPL()

MPLfile.fromHDF(filepath[0])

altrange = np.arange(1000,10000,10)

MPLfile.alt_resample(altrange)

copol = MPLfile.data[0]
crosspol = MPLfile.data[1]

copolvals = np.hstack(copol.values).astype('float32')
crosspolvals = np.hstack(crosspol.values).astype('float32')

depolMPL = crosspolvals/copolvals

depolvals = depolMPL/(depolMPL+1)

copol_mean = np.mean(copolvals)
copol_std = np.std(copolvals)

copol_min = copol_mean-copol_std
copol_max = copol_mean+copol_std

depolhist=h2d.fullhist(depolvals,200,0.24,0.42,-9999.,-8888.)
copolhist=h2d.fullhist(copolvals,200,0.,1.6e-3,-9999.,-8888.)
theOut=h2d.hist2D(copolhist['fullbins'],depolhist['fullbins'],copolhist['numBins'],\
                  depolhist['numBins'])
counts=theOut['coverage']
counts[counts < 1] = 1
logcounts=np.log10(counts)
s = np.shape(counts)
                  
try:
    os.chdir('Plots')
except WindowsError:
    os.makedirs('Plots')
    os.chdir('Plots')

fignum = 0

#fignum+=1
#fig=plt.figure(fignum)
#fig.clf()
#the_axis=fig.add_subplot(111)
#the_axis.plot(depolvals,copolvals,'b+')
#the_axis.set_xlabel('depolvals')
#the_axis.set_ylabel('copolvals')
#the_axis.set_title('raw scatterplot')
#fig.savefig('plot1.png')
#fig.canvas.draw()

cmap=cm.bone
cmap.set_over('r')
cmap.set_under('b')

fignum+=1
fig=plt.figure(fignum)
fig.clf()
axis=fig.add_subplot(111)

im=axis.pcolormesh(depolhist['centers'],copolhist['centers'],logcounts, cmap = cmap)
cb=plt.colorbar(im,extend='both')
title="2-d histogram"
colorbar="log10(counts)"
the_label=cb.ax.set_ylabel(colorbar,rotation=270)
axis.set_xlabel('depolvals')
axis.set_ylabel('copolvals')
axis.set_title(title)
fig.canvas.draw()
fig.savefig('{0}-{1}_hist2d.png'.format(copol.index[0],copol.index[-1]))

fignum+=1
fig = plt.figure(fignum)
fig.clf()
axis = fig.add_subplot(111)
hist = axis.hist(depolvals, bins = 100)

plt.show()

os.chdir(olddir)