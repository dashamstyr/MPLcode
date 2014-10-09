# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 16:23:18 2014

@author: dashamstyr
"""
import random
import numpy as np
import pandas as pan
from matplotlib import pyplot as plt

def randgen(n=100,dist='Gaussian',**kwargs):
    """
    Function to generate an array of random numbers fitting a given PDF
    
    inputs:
    n - number of values to generate
    dist - type of PDF  possible values are:
            'Gaussian': kwargs=mu,sigma
            'Gamma':  kwargs=alpha,beta
            'Lognormal':  kwargs=mu,sigma
            'Uniform': kwargs=minval,maxval
            'Custom': user-defined PDF kwargs=PDF,PDFinputs
    
    outputs:
    randout - an array of length n of random numbers
    """
    mu=kwargs.get('mu',0)
    sigma=kwargs.get('sigma',1)
    alpha=kwargs.get('alpha',1)
    beta=kwargs.get('beta',1)
    minval=kwargs.get('minval',0)
    maxval=kwargs.get('maxval',1)
    PDF=kwargs.get('PDF',[])
    PDFinputs=kwargs.get('PDFinputs',[])
    
    if dist=='Gaussian':
        def tempfun():  return random.gauss(mu,sigma)
    elif dist=='Gamma':
        def tempfun():  return random.gammavariate(alpha,beta)
    elif dist=='Lognormal':
        def tempfun():  return random.lognormvariate(mu,sigma)
    elif dist=='Uniform':
        def tempfun():  return random.uniform(minval,maxval)
    elif dist=='Custom':
        def tempfun(): return PDF(PDFinputs)
    
    randnums=np.empty(n)
    for i in range(n):
        randnums[i]=tempfun()
    
    randout=pan.Series(randnums)
    return randout

def PDFplot(xvals,numbins=25,bounds=[],logplot=False):
    
    if not bounds:
        bounds=(min(xvals),max(xvals))
    binstep=(bounds[1]-bounds[0])/numbins
    binedges=np.arange(bounds[0],bounds[1],binstep)
    print len(bins)
    hvals,binedges=np.histogram(xvals,bins)
    print len(hvals)
    fig1=plt.figure()
    ax1=fig1.add_subplot(111)
    ax1.bar(bins,hvals,width=np.diff(binedges),log=logplot)
    if logplot:
        ax1.set_xscale('log')
    fig1.canvas.draw()

#def PDF 

if __name__=="__main__":
    n=100000
    r=randgen(n=n,dist='Lognormal',mu=500,sigma=1)            
    PDFplot(r.values,numbins=500,logplot=True)    