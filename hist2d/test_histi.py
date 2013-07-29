from __future__ import division
import matplotlib as mpl
import sys,os
import matplotlib.pyplot as plt
from  matplotlib.colorbar import make_axes
import hist2d as h2d
import numpy.random as nr
import numpy as np
import numpy.ma as ma
import imp

invalues={'meanx':450.,
          'stdx':50,
          'meany':-180,
          'stdy':40,
          'rho':0.8}


def makeRandom(meanx=None,stdx=None,meany=None,stdy=None,rho=None):
    nr.seed(50)
    sigma=np.array([stdx**2., rho*stdx*stdy, rho*stdx*stdy, stdy**2.])
    sigma.shape=[2,2]
    meanvec=[meanx,meany]
    outRandom=nr.multivariate_normal(meanvec,sigma,[4000,])
    chan1=outRandom[:,0]
    chan2=outRandom[:,1]
    return (chan1,chan2)

chan1,chan2=makeRandom(**invalues)

bullseye={'meanx':50.,
          'stdx':14,
          'meany':-80,
          'stdy':14,
          'rho':0.0}

chan1B,chan2B=makeRandom(**bullseye)
chan1=np.concatenate((chan1,chan1B))
chan2=np.concatenate((chan2,chan2B))

num1=100
chan1Hist=h2d.fullhist(chan1,num1,0.,650.,-9999,-8888)
num2=75
chan2Hist=h2d.fullhist(chan2,num2,-350.,0.,-9999,-8888)
#put channel 2 on the y axis (rows), so it goes in position 1
theOut=h2d.hist2D(chan2Hist['fullbins'],chan1Hist['fullbins'],chan2Hist['numBins'],\
                  chan1Hist['numBins'])
binList=theOut['indexList']
chan1List=h2d.takeValues(chan1,binList)
chan2List=h2d.takeValues(chan2,binList)
print "chan 1: "

debug=False
if debug:
    for bin,theVec in enumerate(binList):
        row,col=divmod(bin,num2)
        print "index: %d, row %d, column %d" % (bin,row,col)
        print "edges1: ",chan1Hist['edges'][row],chan1Hist['edges'][row+1]
        print "edges2: ",chan2Hist['edges'][col],chan2Hist['edges'][col+1]
        out=chan1List[bin]
        out.sort()
        print out
        out=chan2List[bin]
        out.sort()
        print out
        print '*'*60


fignum=1

fig=plt.figure(fignum)
fig.clf()
the_axis=fig.add_subplot(111)
#for a regular plot, chan2 is the y variable
the_axis.plot(chan1,chan2,'b+')
the_axis.set_xlabel('mean chan 1= %5.0f' % invalues['meanx'])
the_axis.set_ylabel('mean chan 2= %5.0f' % invalues['meany'])
title='raw_coverage'
the_axis.set_title(title)

fignum+=1

fig=plt.figure(fignum)
fig.clf()
the_axis=fig.add_subplot(111)
print theOut['coverage'].shape
#for a countour plot chan2 is the x (row) variable
coverage=theOut['coverage']
#
# suppress Infs when taking the log
#
hit=theOut['coverage'] < 1
coverage[hit]=1.
logcoverage=np.log10(coverage)
coverage=ma.masked_where(logcoverage < 1.e-2,logcoverage)
the_levels=the_axis.contourf(chan1Hist['centers'],chan2Hist['centers'],logcoverage)
axc,kw=make_axes(the_axis,pad=0.02,shrink=0.9,orientation='vertical')
fig.colorbar(the_levels,cax=axc)
the_label=axc.set_ylabel('log10(counts)',rotation=-90)
title='contourf log(coverage)'
the_axis.set_title(title)
fig.canvas.draw()

fignum+=1
fig=plt.figure(fignum)
fig.clf()
the_axis=fig.add_subplot(111)
options={'fignum':fignum}
the_levels=the_axis.imshow(logcoverage)
axc,kw=make_axes(the_axis,pad=0.02,shrink=0.8,orientation='vertical')
fig.colorbar(the_levels,cax=axc)
the_label=axc.set_ylabel('log10(counts)',rotation=-90)
title='imshow log(coverage)'
the_axis.set_title(title)
fig.canvas.draw()

plt.show()
