#cython: embedsignature=True
import numpy as np
from numpy import ma


def getVersion():
    __version__="0.2"
    return "version: %s" % __version__
    

def fullhist(dataVecPy, numbins, mindata, maxdata, missingLowValue, missingHighValue):
    from math import isnan
    """
       given a list or numpy array dataVecPy, bin values in
       numbins between mindata and maxdata, returning a
       python dictionary with edges, centers, counts and
       "fullbins", which is a vector of the same length as
       dataVecPy with the bin of every datapoint
    """
    dataVec=np.ascontiguousarray(dataVecPy,dtype=np.float64)
#    dataPtr=<double*>dataVec.data
    binsize=float(maxdata-mindata)/numbins
    numPts=dataVec.shape[0]
    outcounts = np.zeros([numbins,],dtype=np.int64)
#    countPtr=<Py_ssize_t*> outcounts.data
    bincenters=np.zeros([numbins,],dtype=np.float32)
#    centerPtr=<float*> bincenters.data
    binedges=np.zeros([numbins+1,],dtype=np.float32)
#    edgePtr=<float*> binedges.data
    savebins=np.zeros([numPts,],dtype=np.int64)
#    savebinsPtr=<Py_ssize_t*> savebins.data
    lowcount=0
    highcount=0
    
    for i in range(numPts):
        if isnan(dataVec[i]):
            lowcount+=1
            savebins[i]=missingLowValue
            continue
        else:
            fbin =  int((dataVec[i] - mindata) / binsize)
        if fbin < 0:
            lowcount+=1
            savebins[i]=missingLowValue
            continue
        if fbin > (numbins - 1):
            highcount += 1
            savebins[i] = missingHighValue
            continue
        ibin=fbin
        outcounts[ibin]+=1
        savebins[i]=ibin

    for i in range(numbins + 1):
      binedges[i] = mindata + (i*binsize)
    for i in range(numbins):
      bincenters[i] = float(binedges[i] + binedges[i+1])/2.
    retval={}
    retval["missingLowValue"] = missingLowValue
    retval["missingHighValue"] = missingHighValue
    retval["numBins"] = numbins
    retval["edges"] = binedges
    retval["centers"] = bincenters
    retval["counts"] = outcounts
    retval["lowcounts"] = lowcount
    retval["highcounts"] = highcount
    retval["fullbins"] = savebins
    return retval


def hist2D(xBinPy,yBinPy,numXbins,numYbins):
    """
      xBinArray is a vector of bin indices, each pixel gets a bin number
      yBinArray is a vector of bin indices, each pixel gets a bin number
      numXbins is the total number of bin indices for x
      numYbins is the total number of bin indices for y
      converageMap is a 2-d histogram with the number of points
      in each 2d bin
    """
    xBinArray=np.ascontiguousarray(xBinPy,dtype=np.int64)
    yBinArray=np.ascontiguousarray(yBinPy,dtype=np.int64)
    numXDataPoints=xBinArray.shape[0]
    numYDataPoints=yBinArray.shape[0]
    if numXDataPoints != numYDataPoints:
        raise ValueError('need x and y fields of equal size')
    numBins2D=numXbins*numYbins
    binVecs=[]
    for i in range(numBins2D):
      binVecs.append([])
    x=0
    y=0
    index=0
    #drop the indexes into a nensted list
    for i in range(numXDataPoints):
      if (xBinArray[i] > -1) & (yBinArray[i] > -1):
        x = xBinArray[i]
        y = yBinArray[i]
        #2D row major, numYbins is number of columns, numXbins is number of rows
        #if numXbins=10 and numYbins=5, then an (x,y) of (5,3) gives
        #an index of 28
        index = numYbins*x + y
        binVecs[index].append(i)
    #return an 2D numpy array with the number of points in each cell
    coverageMap = np.zeros(numBins2D,dtype=np.int64)
    #convert list of list to list of np.arrays
    arrayList=[]
    for i in range(numBins2D):
      arrayList.append(np.array(binVecs[i],dtype=np.int64))
      #number of pixels in each bin
      coverageMap[i]=len(binVecs[i])
     
    coverageMap = np.reshape(coverageMap,(numXbins,numYbins))
    retval={}
    retval["coverage"]= coverageMap
    retval["indexList"]= arrayList
    return retval


def takeValues(dataVectorPy, indexList):
    """
    do a take of the indices in indexList to populate a new list of data
    filled with dataVector values.  See findMean and findMedian
    below for usage
    """
    dataVector=np.ascontiguousarray(dataVectorPy,dtype=np.float32)
    dataVector=dataVector.reshape(-1)
    outList=[]
    numBins2D=len(indexList)
    for i in range(numBins2D):
        indexVec=indexList[i]
        numDataPoints=len(indexVec)
        takeVec=np.zeros([numDataPoints,],dtype=np.float32)
        for j in range(numDataPoints):
          takeVec[j]=dataVector[indexVec[j]]
        outList.append(takeVec)
    return outList

def findMean(dataVectorPy,indexList,maskedValue= -9999.):
    """
    find the mean of binned variables
    """
    dataVector=np.ascontiguousarray(dataVectorPy,dtype=np.float32)
    dataVector=dataVector.reshape(-1)
    dataList=takeValues(dataVector,indexList)
    dataCount=len(dataList)
    outList=[]
    areaWeightedOutList=[]
    gridCounts=[]
    for i in range(dataCount):
        theData=dataList[i]
        if len(theData) > 0:
            outList.append(theData.mean())
            areaWeightedOutList.append(theData.mean()*len(theData))
            #print "appending: ",len(theData)
            gridCounts.append(len(theData))
        else:
            outList.append(maskedValue)
            areaWeightedOutList.append(maskedValue)
            gridCounts.append(0)
    outVec=np.array(outList,dtype=np.float32)
    outAreaWeightedVec=np.array(areaWeightedOutList,dtype=np.float32)
    outArray=ma.masked_where(outVec==maskedValue,outVec)    
    outAreaArray=ma.masked_where(outAreaWeightedVec==maskedValue,outAreaWeightedVec)
    gridCounts=np.array(gridCounts,dtype=np.int64)
    return (outArray,outAreaArray,gridCounts)

def findMedian(dataVector,indexList,maskedValue= -9999.):
    """
    do a take of the indices in indexList to populate a new list of data
    filled with dataVector values
    """
    dataList=takeValues(dataVector,indexList)
    dataCount=len(dataList)
    outList=[]
    for i in range(dataCount):
        theData=dataList[i]
        if len(theData) > 0:
            outList.append(np.median(theData))
        else:
            outList.append(maskedValue)
    outVec=np.array(outList,dtype=np.float32)
    outArray=ma.masked_where(outVec==maskedValue,outVec)    
    return outArray

def althist(datavalsPy, altvalsPy, numdatbins, histrange=None):
    """
      datavalsPy is an array of data values with each column representing an altitude
      altvalsPy is a vector of altitudes with one value for each column
      valspecol is the number of data values per column
      histrange = (min,max) a tuple representing the minimim and maximum data values
      converageMap is a 2-d histogram with the number of points
      in each 2d bin
    """
    datavals=np.ascontiguousarray(datavalsPy)
    altvals=np.ascontiguousarray(altvalsPy)
    
    if numdatbins >= datavals.shape[0]:
        raise ValueError("Number of bins cannot exceed data set")
    else:
        numXDataPoints = numdatbins 
    
    numYDataPoints=altvals.shape[0]
    if datavals.shape[1] != altvals.shape[0]:
        raise ValueError('need one altitude per data column')
#    numBins2D=numXDataPoints*numYDataPoints
    
    coverageMap = np.empty([numXDataPoints,numYDataPoints])
    indexArray = np.empty([numXDataPoints+1,numYDataPoints])
    #drop the indexes into a nested list
    for i in range(numYDataPoints):
        coverageMap[:,i],indexArray[:,i] = np.histogram(datavals[:,i],bins=numdatbins,range=histrange)
    #      if (xBinArray[i] > -1) & (yBinArray[i] > -1):
    #        x = xBinArray[i]
    #        y = yBinArray[i]
    #        #2D row major, numYbins is number of columns, numXbins is number of rows
    #        #if numXbins=10 and numYbins=5, then an (x,y) of (5,3) gives
    #        #an index of 28
    #        index = numYbins*x + y
    #        binVecs[index].append(i)
    #return an 2D numpy array with the number of points in each cell
    #    coverageMap = np.zeros(numBins2D,dtype=np.int64)
    #    #convert list of list to list of np.arrays
    #    arrayList=[]
    #    for i in range(numBins2D):
    #      arrayList.append(np.array(binVecs[i],dtype=np.int64))
    #      #number of pixels in each bin
    #      coverageMap[i]=len(binVecs[i])
    #     
    #    coverageMap = np.reshape(coverageMap,(numXbins,numYbins))
        
    retval={}
    retval["coverage"]= coverageMap
    retval["indexList"]= indexArray
    return retval