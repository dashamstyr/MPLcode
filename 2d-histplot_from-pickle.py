# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:10:30 2013

@author: dashamstyr
"""

import site
site.addsitedir('/Users/phil/repos/MPLcode')
site.addsitedir('/Users/phil/lib/python')
import MPLtools as mtools
import os
import numpy as np
import pandas as pan
import matplotlib.pyplot as plt
import pickle
from matplotlib.colors import Normalize
from matplotlib import cm

#from plot_fields import loadvars

## import fasthist as fh
## from binit import fastin

from matplotlib.colors import Normalize
from matplotlib import cm

#os.chdir('C:\Program Files (x86)\SigmaMPL\DATA')
#
#filepath = mtools.get_files('Select MPL file', filetype = ('.h5', '*.h5'))
#
#MPLfile = mtools.MPL()
#
#MPLfile.fromHDF(filepath[0])
#
#altrange = np.arange(150,4000,10)
#
#MPLfile.alt_resample(altrange)

#copol = MPLfile.data[0]
#crosspol = MPLfile.data[1]

MPLdat = pickle.load(open('201304222100-201304222300_proc.pickle','rb'))

copol = MPLdat[0]
crosspol = MPLdat[1]

copolvals = np.hstack(copol.values).astype('float32')
crosspolvals = np.hstack(crosspol.values).astype('float32')

depolMPL = crosspolvals/copolvals

depolvals = depolMPL/(depolMPL+1)

copol_mean = np.mean(copolvals)
copol_std = np.std(copolvals)

copol_min = copol_mean-copol_std
copol_max = copol_mean+copol_std

fignum=1
fig=plt.figure(fignum)
fig.clf()
the_axis=fig.add_subplot(111)
the_axis.plot(depolvals,copolvals,'b+')
the_axis.set_xlabel('depolvals')
the_axis.set_ylabel('copolvals')
the_axis.set_title('raw scatterplot')
fig.savefig('plot1.png')
fig.canvas.draw()

import fasthist as h2d
depolhist=h2d.fullhist(depolvals,20,0.24,0.42,-9999.,-8888.)
copolhist=h2d.fullhist(copolvals,20,0.,1.6e-3,-9999.,-8888.)
theOut=h2d.hist2D(copolhist['fullbins'],depolhist['fullbins'],copolhist['numBins'],\
                  depolhist['numBins'])

cmap=cm.bone
cmap.set_over('r')
cmap.set_under('b')

fig=plt.figure(2)
fig.clf()
axis1=fig.add_subplot(111)
counts=theOut['coverage']
counts[counts < 1] = 1
logcounts=np.log10(counts)
im=axis1.pcolormesh(depolhist['centers'],copolhist['centers'],logcounts,
                    cmap=cmap)
cb=plt.colorbar(im,extend='both')
title="2-d histogram"
colorbar="log10(counts)"
the_label=cb.ax.set_ylabel(colorbar,rotation=270)
axis1.set_xlabel('depolvals')
axis1.set_ylabel('copolvals')
axis1.set_title(title)
fig.canvas.draw()
fig.savefig('plot2.png')
plt.show()

