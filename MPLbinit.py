#updated 2014/May/22
from __future__ import division
#see http://pysclint.sourceforge.net/pyhdf/
import numpy as np
import matplotlib.pyplot as plt
import glob, sys
from types import MethodType

def progress(datanum,tot_loops):
    """
       print a message about percent complete
       input:      datanum: current index (int)
       tot_loops:  index when loop finishes
    """   
    the_frac=np.int(datanum/tot_loops*100.)
    sys.stdout.write("\rpercent complete: %d%%" % the_frac)
    sys.stdout.flush()


class binit(object):
    """histograms a vector of data, returning the bin index of every
       datapoint

       Constructor
       -----------

          bin_obj=binit(minval,maxval,numbins,
                         missingLowValue,missingHighValue)

        Parameters
        ----------

         minval: float
            left edge of smallest bin
         maxval: float
            right edge of largest bin
         numbins: int
            number of bins
         missingLowValue: float
            bin number indicating data is smaller than minval
         missingHighValue: float
            bin number indicating data is larger than maxval
                         
    """
    def __init__(self,minval,maxval,numbins,missingLowValue=-99999,missingHighValue=99999):
        
        """
        inputs:
        missingLowValue = integer substituted for values below minval
        missingHighValue = integer substituted for values above maxval
        minval = minimum value
        maxval = maximum value
        numbins = number of bins to break vector into
        binsize = size of individual bins
        minval,maxval and either numbins or binsize must be defined
        """

        self.minval=minval
        self.maxval=maxval
        self.numbins=numbins 
        self.missingLowValue=missingLowValue
        self.missingHighValue=missingHighValue        
        self.bin_edges=np.arange(self.minval,self.maxval + 0.5*self.binsize,self.binsize)            
        self.bin_centers= (self.bin_edges[:-1] + self.bin_edges[1:])/2.
        self.binsize = (maxval-minval)/numbins
        self.data=[]
        
                    
    def do_bins(self,data_vec):
        """
           bin the items in data_vec into self.numbins
           see binit docstring for details

           parameters
           ----------

              data_vec: numpy 1d array (float)
                 vector of floating point numbers to be binned

           returns:
           --------
             bin_count: numpy vector of length self.numbins (int)
                vector containing bin counts
             
             bin_index: numpy vector of len(data_vec)
                vector containing the bin number of each datapoint

             lowcount: int
                number of points that were smaller
                than the smallest bin
                
             highcount: int
                number of points that were larger than the largest
                bin

           example
           -------
             
             (bin_count,bin_index,lowcount,hightcount)=obj.do_bins(lat_vec) 
        """   
        bin_index=np.empty_like(data_vec,dtype=np.int)
        bin_count=np.zeros([self.numbins],dtype=np.int)
        lowcount=0
        highcount=0
        tot_loops=len(data_vec)
        for datanum,dataval in enumerate(data_vec):
            if np.mod(datanum,10000)==0:
                progress(datanum,tot_loops)
            float_bin =  ((dataval - self.minval) /self.binsize)
            if float_bin < 0:
                lowcount+=1
                bin_index[datanum]=self.missingLowValue
                continue
            if float_bin > self.numbins:
                highcount += 1
                bin_index[datanum] = self.missingHighValue
                continue
            ibin=int(float_bin)
            bin_count[ibin]+=1
            bin_index[datanum]=ibin
        
        self.bin_count = bin_count
        self.bin_index = bin_index
        self.lowcount = lowcount
        self.highcount = highcount
        self.data = data_vec
        
        return(self.bin_count,self.bin_index,self.lowcount,self.highcount)

    def get_centers(self):
        """
          Get the bin centers for the historgram
        """
        return self.bin_centers

    def get_edges(self):
        """
          Get the bin edges for the historgram
        """
        return self.bin_edges
        
    def append(self,newbinnit,method='union',average=False):
        """
        adds a new binnit object to the current one according to the prescribed 
        method:
        
        left = original maxval,minval,binsize used
        right = new maxval,minval,binsize used
        union = union of both bins used - smaller of two binsizes used
        intersection = intersection of both used - smaller of two binsizes used
        
        average = if ture, values are averaged together instead of summed
        
        """
        if method=="right":
            self.minval=newbinnit.minval
            self.maxval=newbinnit.maxval
            self.numbins=newbinnit.numbins
            self.binsize=newbinnit.binsize
            self.bin_edges=newbinnit.bin_edges
            self.bin_centers= newbinnit.bin_centers
        elif method=="union":
            self.minval=np.min([self.minval,newbinnit.minval])
            self.maxval=np.max([self.maxval,newbinnit.maxval])
            self.binsize=np.min([self.binsize,newbinnit.binsize])
            self.bin_edges=np.arange(self.minval,self.maxval + 0.5*self.binsize,self.binsize)            
            self.bin_centers= (self.bin_edges[:-1] + self.bin_edges[1:])/2.
            self.numbins=int((self.maxval-self.minval)/self.binsize)
        elif method=="intersection":
            self.minval=np.max([self.minval,newbinnit.minval])
            self.maxval=np.min([self.maxval,newbinnit.maxval])
            self.binsize=np.min([self.binsize,newbinnit.binsize])
            self.bin_edges=np.arange(self.minval,self.maxval + 0.5*self.binsize,self.binsize)            
            self.bin_centers= (self.bin_edges[:-1] + self.bin_edges[1:])/2.
            self.numbins=int((self.maxval-self.minval)/self.binsize)
                
        data_vec = self.data+newbinnit.data                    
        bin_count,bin_index,lowcount,highcount = self.do_bins(data_vec)
#
#class fastbin(binit):
#    def __init__(self,minval,maxval,numbins,missingLowValue,missingHighValue):
#        binit.__init__(self,minval,maxval,numbins,missingLowValue,missingHighValue)
#
##
## replace do_bins with the cython version from fastbinit.pyx
##
#from fastbinit import do_bins
#fastbin.do_bins=MethodType(do_bins, None, fastbin)

    



        
