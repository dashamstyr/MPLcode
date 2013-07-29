#cython: embedsignature=True
import numpy as np
cimport numpy as np

from numpy import ma


def getVersion():
    __version__="0.2"
    return "version: %s" % __version__
    

def  fullhist(object dataVecPy, int numbins, float mindata, 
              float maxdata,int missingLowValue,int missingHighValue):
    """
       given a list or numpy array dataVecPy, bin values in
       numbins between mindata and maxsdata, returning a
       python dictionary with edges, centers, counts and
       "fullbins", which is a vector of the same length as
       dataVecPy with the bin of every datapoint
    """
    dataVecPy=np.ascontiguousarray(dataVecPy,dtype=np.float64)
    cdef np.ndarray[np.double_t,ndim=1] dataVec = dataVecPy
    cdef double* dataPtr=<double*>dataVec.data
    cdef float binsize=(maxdata-mindata)/numbins
    cdef np.int64_t numPts=dataVec.shape[0]
    cdef np.ndarray[np.int64_t,ndim=1] outcounts = np.zeros([numbins,],dtype=np.int64)
    cdef Py_ssize_t* countPtr=<Py_ssize_t*> outcounts.data
    cdef np.ndarray[np.float32_t,ndim=1] bincenters=np.zeros([numbins,],dtype=np.float32)
    cdef float* centerPtr=<float*> bincenters.data
    cdef np.ndarray[np.float32_t,ndim=1] binedges=np.zeros([numbins+1,],dtype=np.float32)
    cdef float* edgePtr=<float*> binedges.data
    cdef np.ndarray[np.int64_t,ndim=1] savebins=np.zeros([numPts,],dtype=np.int64)
    cdef Py_ssize_t* savebinsPtr=<Py_ssize_t*> savebins.data
    cdef float fbin
    cdef long lowcount=0, highcount=0
    cdef int i,ibin
    for i in range(numPts):
        fbin =  ((dataPtr[i] - mindata) / binsize)
        if fbin < 0:
            lowcount+=1
            savebinsPtr[i]=missingLowValue
            continue
        if fbin > (numbins - 1):
            highcount += 1
            savebinsPtr[i] = missingHighValue
            continue
        ibin=<int>fbin
        countPtr[ibin]+=1
        savebinsPtr[i]=ibin

    for i in range(numbins + 1):
      edgePtr[i] = mindata + (i*binsize)
    for i in range(numbins):
      centerPtr[i] = (edgePtr[i] + edgePtr[i+1])/2.
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


def hist2D(object xBinPy,object yBinPy,int numXbins, int numYbins):
    """
      xBinArray is a vector of bin indices, each pixel gets a bin number
      yBinArray is a vector of bin indices, each pixel gets a bin number
      numXbins is the total number of bin indices for x
      numYbins is the total number of bin indices for y
      converageMap is a 2-d histogram with the number of points
      in each 2d bin
    """
    xBinPy=np.ascontiguousarray(xBinPy,dtype=np.int64)
    cdef np.ndarray[np.int64_t,ndim=1] xBinArray = xBinPy
    yBinPy=np.ascontiguousarray(yBinPy,dtype=np.int64)
    cdef np.ndarray[np.int64_t,ndim=1] yBinArray = yBinPy
    cdef Py_ssize_t* xdataPtr=<Py_ssize_t*>xBinArray.data
    cdef Py_ssize_t* ydataPtr=<Py_ssize_t*>yBinArray.data
    cdef np.int64_t  numXDataPoints=xBinArray.shape[0]
    cdef np.int64_t numYDataPoints=yBinArray.shape[0]
    if numXDataPoints != numYDataPoints:
        raise ValueError('need x and y fields of equal size')
    cdef int numBins2D=numXbins*numYbins
    binVecs=[]
    cdef int i
    for i in range(numBins2D):
      binVecs.append([])
    cdef int x = 0, y=0, index=0
    #drop the indexes into a nensted list
    for i in range(numXDataPoints):
      if (xdataPtr[i] > -1) & (ydataPtr[i] > -1):
        x = xdataPtr[i]
        y = ydataPtr[i]
        #2D row major, numYbins is number of columns, numXbins is number of rows
        #if numXbins=10 and numYbins=5, then an (x,y) of (5,3) gives
        #an index of 28
        index = numYbins*x + y
        binVecs[index].append(i)
    #return an 2D numpy array with the number of points in each cell
    cdef np.ndarray[np.int64_t,ndim=2] coverageMap = np.zeros([numXbins,numYbins],dtype=np.int64)
    cdef Py_ssize_t* coveragePtr = <Py_ssize_t*> coverageMap.data
    cdef int numPoints
    #convert list of list to list of np.arrays
    arrayList=[]
    for i in range(numBins2D):
      arrayList.append(np.array(binVecs[i],dtype=np.int64))
      #number of pixels in each bin
      coveragePtr[i]=len(binVecs[i])
    retval={}
    retval["coverage"]= coverageMap
    retval["indexList"]= arrayList
    return retval


def takeValues(object dataVectorPy, object indexList):
    """
    do a take of the indices in indexList to populate a new list of data
    filled with dataVector values.  See findMean and findMedian
    below for usage
    """
    dataVectorPy=np.ascontiguousarray(dataVectorPy,dtype=np.float32)
    dataVectorPy=dataVectorPy.reshape(-1)
    cdef np.ndarray[np.float32_t,ndim=1] dataVector = dataVectorPy
    outList=[]
    cdef float* dataPtr=<float*>dataVector.data
    cdef np.ndarray[np.int64_t,ndim=1] indexVec
    cdef np.ndarray[np.float32_t,ndim=1] takeVec
    cdef int numBins2D=len(indexList)
    cdef int numDataPoints,i,j
    cdef float* takePtr
    for i in range(numBins2D):
        indexVec=indexList[i]
        numDataPoints=len(indexVec)
        takeVec=np.zeros([numDataPoints,],dtype=np.float32)
        takePtr=<float*>takeVec.data
        for j in range(numDataPoints):
          takePtr[j]=dataPtr[indexVec[j]]
        outList.append(takeVec)
    return outList

def findMean(object dataVectorPy, object indexList,maskedValue= -9999.):
    """
    find the mean of binned variables
    """
    dataVectorPy=np.ascontiguousarray(dataVectorPy,dtype=np.float32)
    dataVectorPy=dataVectorPy.reshape(-1)
    cdef np.ndarray[np.float32_t,ndim=1] dataVector = dataVectorPy
    dataList=takeValues(dataVector,indexList)
    cdef int dataCount=len(dataList)
    cdef int i
    outList=[]
    areaWeightedOutList=[]
    gridCounts=[]
    cdef np.ndarray[float,ndim=1] theData
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

def findMedian(np.ndarray[float,ndim=1] dataVector, object indexList,maskedValue= -9999.):
    """
    do a take of the indices in indexList to populate a new list of data
    filled with dataVector values
    """
    dataList=takeValues(dataVector,indexList)
    cdef int dataCount=len(dataList)
    cdef int i
    outList=[]
    cdef np.ndarray[float,ndim=1] theData
    for i in range(dataCount):
        theData=dataList[i]
        if len(theData) > 0:
            outList.append(np.median(theData))
        else:
            outList.append(maskedValue)
    outVec=np.array(outList,dtype=np.float32)
    outArray=ma.masked_where(outVec==maskedValue,outVec)    
    return outArray

## def findTemp(cnp.ndarray theta, cnp.ndarray press):
##     tempout=np.zeros([theta.size,])
##     cdef float* tempPtr=<float*>tempout.data
##     cdef int numRows=theta.shape[0]
##     cdef int numCols=theta.shape[1]
##     theta=theta.ravel()
##     press=press.ravel()
##     cdef int i
##     cdef int numPts=theta.size
##     cdef cnp.ndarray logtheta=np.log(theta)
##     cdef float logP0=np.log(1.e5)
##     cdef cnp.ndarray logP=np.log(press)
##     cdef float logTemp
##     cdef float RdoCp=287./1004.
##     print "in loop findTemp"
##     for i in range(numPts):
##         tempout[i]=logtheta[i] - RdoCp*(logP0 - logP[i])
##     tempout=np.exp(tempout)
##     return tempout
