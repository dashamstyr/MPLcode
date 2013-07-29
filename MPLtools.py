 # -*- coding: utf-8 -*-
"""
mpltools.py
A bag of tools to be used in processing and interpreting MPL class data
collected by miniMPL
Created on Wed Apr 24 12:08:57 2013

@author: Paul Cottle

"""

import numpy as np
import time, datetime
import array
import pandas as pan
from scipy import constants as const
from scipy.interpolate import interp1d
import tables
from Tkinter import Tk
import tkFileDialog
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm

def set_dir(titlestring):
    
    master = Tk()
    master.withdraw() #hiding tkinter window
     
    file_path = tkFileDialog.askdirectory(title=titlestring)
     
    if file_path != "":
       return str(file_path)
     
    else:
       print "you didn't open anything!"
     
def get_files(titlestring,filetype = ('.txt','*.txt')):
         
    master = Tk()
    master.withdraw() #hiding tkinter window

    file_path = []
     
    file_path = tkFileDialog.askopenfilename(title=titlestring, filetypes=[filetype,("All files",".mpl")],multiple='True')
     
    if file_path != "":
       file_path = str(file_path)
       file_path = file_path.lstrip('{')
       file_path = file_path.rstrip('}')
       
       return file_path.split('} {')
    else:
       print "you didn't open anything!"
    
    
def MPLtoHDF(filename, appendflag = 'False'):
    
    #Define class header to create columns to hold header data from .mpl binary file
    class header(IsDescription):
        unitnum = tables.UInt32Col(1)
        version = tables.UInt32Col(1)
        
        timestamp = tables.Time32Col(1)
        
        shotsum = tables.UInt32Col(1)  #total number of shots collected per profile
        trigfreq = tables.UInt32Col(1) #laser trigger frequency (usually 2500 Hz)
        energy = tables.UInt32Col(1) #mean of laser energy monitor in uJ
        temp_0 = tables.UInt32Col(1)  #mean of A/D#0 readings*100
        temp_1 = tables.UInt32Col(1)  #mean of A/D#1 readings*100
        temp_2 = tables.UInt32Col(1)  #mean of A/D#2 readings*100
        temp_3 = tables.UInt32Col(1)  #mean of A/D#3 readings*100
        temp_4 = tables.UInt32Col(1)  #mean of A/D#4 readings*100
        
        bg_avg1 = tables.Float32Col(1) #mean background signal value for channel 1
        bg_std1 = tables.Float32Col(1) #standard deviation of backgruond signal for channel 1

        numchans = tables.UInt16Col(1) #number of channels
        numbins = tables.UInt32Col(1) #total number of bins per channel

        bintime = tables.Float32Col(1)  #bin width in seconds
        
        rangecal = tables.Float32Col(1) #range offset in meters, default is 0

        databins = tables.UInt16Col(1)  #number of bins not including those used for background
        scanflag = tables.UInt16Col(1)  #0: no scanner, 1: scanner
        backbins = tables.UInt16Col(1)  #number of background bins

        az = tables.Float32Col(1)  #scanner azimuth angle
        el = tables.Float32Col(1)  #scanner elevation angle
        deg = tables.Float32Col(1) #compass degrees (currently unused)
        pvolt0 = tables.Float32Col(1) #currently unused
        pvolt1 = tables.Float32Col(1) #currently unused
        gpslat = tables.Float32Col(1) #GPS latitude in decimal degreees (-999.0 if no GPS)
        gpslon = tables.Float32Col(1) #GPS longitude in decimal degrees (-999.0 if no GPS)
        cloudbase = tables.Float32Col(1) #cloud base height in [m]

        baddat = tables.BoolCol(1)  #0: good data, 1: bad data
        version = tables.BoolCol(1) #version of file format.  current version is 1

        bg_avg2 = tables.Float32Col(1) #mean background signal for channel 2
        bg_std2 = tables.Float32Col(1) #mean background standard deviation for channel 2

        mcs = tables.BoolCol(8)  #MCS mode register  Bit#7: 0-normal, 1-polarization
                                     #Bit#6-5: polarization toggling: 00-linear polarizer control
                                     #01-toggling pol control, 10-toggling pol control 11-circular pol control

        firstbin = tables.UInt16Col(1)  #bin # of first return data
        systype = tables.BoolCol(1)   #0: standard MPL, 1: mini MPL
        syncrate = tables.UInt16Col(1)  #mini-MPL only, sync pulses seen per second
        firstback = tables.UInt16Col(1) #mini-MPL only, first bin used for background calcs
        headersize2 = tables.UInt16Col(1) #size of additional header data (currently unused)
    
    
    
    h5filename = filename.split('.')[0]+'.h5'
    with tables.openFile(h5filename, mode = 'w', title = 'MPL data file') as h5file:
        
        headertbl = h5file.create_table('/','Header',header,'Ancillary Data')
          
        headerdat = headertbl.row
        
        binfile = open(filename,'rb')
            
        intarray16 = array.array('H')
        intarray32 = array.array('L')
        floatarray = array.array('f')
        byte_array = array.array('B')
        datavals = array.array('f')
        
        
        intarray16.fromfile(binfile, 8)
        intarray32.fromfile(binfile, 8)
        floatarray.fromfile(binfile, 2)
        intarray16.fromfile(binfile, 1)
        intarray32.fromfile(binfile, 1)
        floatarray.fromfile(binfile, 2)
        intarray16.fromfile(binfile, 3)
        floatarray.fromfile(binfile, 8)
        byte_array.fromfile(binfile, 2)
        floatarray.fromfile(binfile, 2)
        byte_array.fromfile(binfile, 1)
        intarray16.fromfile(binfile, 1)
        byte_array.fromfile(binfile, 1)
        intarray16.fromfile(binfile, 3)
        
        headerdat['unitnum'] = intarray16[0]
        headerdat['version'] = intarray16[1]
        year = intarray16[2]
        month = intarray16[3]
        day = intarray16[4]
        hour = intarray16[5]
        minute = intarray16[6]
        second = intarray16[7]
        dt = datetime.datetime(year,month,day,hour,minute,second)
    
        headerdat['timestamp'] = time.mktime(dt.timetuple())
    
        headerdat['shotsum'] = intarray32[0]  #total number of shots collected per profile
        headerdat['trigfreq'] = intarray32[1] #laser trigger frequency (usually 2500 Hz)
        headerdat['energy'] = intarray32[2]/1000.0  #mean of laser energy monitor in uJ
        headerdat['temp_0'] = intarray32[3]/1000.0  #mean of A/D#0 readings*100
        headerdat['temp_1'] = intarray32[4]/1000.0  #mean of A/D#1 readings*100
        headerdat['temp_2'] = intarray32[5]/1000.0  #mean of A/D#2 readings*100
        headerdat['temp_3'] = intarray32[6]/1000.0  #mean of A/D#3 readings*100
        headerdat['temp_4'] = intarray32[7]/1000.0  #mean of A/D#4 readings*100
        
        headerdat['bg_avg1'] = floatarray[0] #mean background signal value for channel 1
        headerdat['bg_std1'] = floatarray[1] #standard deviation of backgruond signal for channel 1
    
        headerdat['numchans'] = intarray16[8] #number of channels
        headerdat['numbins'] = intarray32[8] #total number of bins per channel
    
        headerdat['bintime'] = floatarray[2]  #bin width in seconds
        
        headerdat['rangecal'] = floatarray[3] #range offset in meters, default is 0
    
        headerdat['databins'] = intarray16[9]  #number of bins not including those used for background
        headerdat['scanflag'] = intarray16[10]  #0: no scanner, 1: scanner
        headerdat['backbins'] = intarray16[11]  #number of background bins
    
        headerdat['az'] = floatarray[4]  #scanner azimuth angle
        headerdat['el'] = floatarray[5]  #scanner elevation angle
        headerdat['deg'] = floatarray[6] #compass degrees (currently unused)
        headerdat['pvolt0'] = floatarray[7] #currently unused
        headerdat['pvolt1'] = floatarray[8] #currently unused
        headerdat['gpslat'] = floatarray[9] #GPS latitude in decimal degreees (-999.0 if no GPS)
        headerdat['gpslon'] = floatarray[10]#GPS longitude in decimal degrees (-999.0 if no GPS)
        headerdat['cloudbase'] = floatarray[11] #cloud base height in [m]
    
        headerdat['baddat'] = byte_array[0]  #0: good data, 1: bad data
        headerdat['version'] = byte_array[1] #version of file format.  current version is 1
    
        headerdat['bg_avg2'] = floatarray[12] #mean background signal for channel 2
        headerdat['bg_std2'] = floatarray[13] #mean background standard deviation for channel 2
    
        headerdat['mcs'] = byte_array[2]  #MCS mode register  Bit#7: 0-normal, 1-polarization
                                     #Bit#6-5: polarization toggling: 00-linear polarizer control
                                     #01-toggling pol control, 10-toggling pol control 11-circular pol control
    
        headerdat['firstbin'] = intarray16[12]  #bin # of first return data
        headerdat['systype'] = byte_array[3]   #0: standard MPL, 1: mini MPL
        headerdat['syncrate'] = intarray16[13]  #mini-MPL only, sync pulses seen per second
        headerdat['firstback'] = intarray16[14] #mini-MPL only, first bin used for background calcs
        headerdat['headersize2'] = intarray16[15] #size of additional header data (currently unused)
                 
        numbins = headerdat['numbins']
        numchans = headerdat['numchans']         
    
        while True:
            try:
                datavals.fromfile(binfile,(numbins*numchans))                    
            except EOFError:
                break
                
            altstep = header['bintime']*const.c #altitude step in meters
            maxalt = numbins*altstep
            minalt = header['rangecal']
            altrange = np.arange(minalt,maxalt,altstep,dtype='float')
            
            timestep = datetime.timedelta(seconds = header['shotsum']/header['trigfreq']) #time between profiles in seconds
            dt = datetime.fromtimestamp(header['timestamp'])
      
    if numchans == 2:
        profdat_copol = []
        profdat_crosspol = []
        
        n = 0
        while n < len(datavals):
            profdat_crosspol.append(datavals[n:n+numbins]) 
            n += numbins                               
            profdat_copol.append(datavals[n:n+numbins])
            
            n += numbins + 32
            
            #note: this more adaptive algorithm should have worked
            #I'm not sure why it doesn't so I'm going with the
            #patch for now - 04-25-13
#                            m = n+numbins
#                            temp_prof_crosspol = datarray[n:m]
#                            temp_prof_copol = datarray[m:m+numbins]
#                            
#                            
#                            if all(temp_prof_crosspol > 1e-6):                             
#                                profdat_crosspol.append(temp_prof_crosspol)                                
#                                profdat_copol.append(temp_prof_copol)
#                                n += 2*numbins
#                            else:
#                                n += 1
        copoldat = np.array(profdat_copol)
        crosspoldat = np.array(profdat_crosspol)
        numsteps = np.shape(copoldat)[0]

        timerange = []
        for t in range(0,numsteps):
            timerange.append(dt)
            dt += timestep  
            
        df_copol = pan.DataFrame(copoldat,index = timerange,columns = altrange)
        df_crosspol = pan.DataFrame(crosspoldat,index = timerange,columns = altrange)
        
        df_copol.to_hdf(h5filename,'copol', append = appendflag)
        df_crosspol.to_hdf(h5filename,'crosspol', append = appendflag)
        
    else:
        raise ValueError('Wrong number of channels')
    
       

class MPL:
    """
    This is a class type generated by unpacking a binary file generated by
    the mini-MPL lidar
    
    It includes two subclasses: header and data
    The metadata in the header is described in the MPL manual pp 37-38
    The data consists of a 2-D array of floats separated into channels
    
    copol = data measured in laser polarization
    crosspol = data measured in cross polarization
    
    """
#    import numpy as np
#    import datetime
#    import array
#    import pandas as pan
#    from scipy import constants as const
#    from scipy.interpolate import interp1d
#    import tables
#    from copy import deepcopy
#    import matplotlib.pyplot as plt
#    from matplotlib.colors import Normalize
#    from matplotlib import cm
        
    def __init__(self,filename=[]):
        
        self.data = [] #slot for lidar data array
        
        class header:
            def __init__(self):
                self.unitnum = []
                self.version = []
                self.datetime = []
                self.shotsum = []  #total number of shots collected per profile
                self.trigfreq = [] #laser trigger frequency (usually 2500 Hz)
                self.energy = [] #mean of laser energy monitor in uJ
                self.temp_0 = []  #mean of A/D#0 readings*100
                self.temp_1 = []  #mean of A/D#1 readings*100
                self.temp_2 = []  #mean of A/D#2 readings*100
                self.temp_3 = []  #mean of A/D#3 readings*100
                self.temp_4 = []  #mean of A/D#4 readings*100
                
                self.bg_avg1 = [] #mean background signal value for channel 1
                self.bg_std1 = [] #standard deviation of backgruond signal for channel 1

                self.numchans = [] #number of channels
                self.numbins = [] #total number of bins per channel

                self.bintime = []  #bin width in seconds
                
                self.rangecal = [] #range offset in meters, default is 0

                self.databins = []  #number of bins not including those used for background
                self.scanflag = []  #0: no scanner, 1: scanner
                self.backbins = []  #number of background bins

                self.az = []  #scanner azimuth angle
                self.el = []  #scanner elevation angle
                self.deg = [] #compass degrees (currently unused)
                self.pvolt0 = [] #currently unused
                self.pvolt1 = [] #currently unused
                self.gpslat = [] #GPS latitude in decimal degreees (-999.0 if no GPS)
                self.gpslon = []#GPS longitude in decimal degrees (-999.0 if no GPS)
                self.cloudbase = [] #cloud base height in [m]

                self.baddat = []  #0: good data, 1: bad data
                self.version = [] #version of file format.  current version is 1

                self.bg_avg2 = [] #mean background signal for channel 2
                self.bg_std2 = [] #mean background standard deviation for channel 2

                self.mcs = []  #MCS mode register  Bit#7: 0-normal, 1-polarization
                                             #Bit#6-5: polarization toggling: 00-linear polarizer control
                                             #01-toggling pol control, 10-toggling pol control 11-circular pol control

                self.firstbin = []  #bin # of first return data
                self.systype = []   #0: standard MPL, 1: mini MPL
                self.syncrate = []  #mini-MPL only, sync pulses seen per second
                self.firstback = [] #mini-MPL only, first bin used for background calcs
                self.headersize2 = [] #size of additional header data (currently unused)
        
        self.header = header()
    
    def fromHDF(self, filename):
        
        with tables.openFile(filename,'r+') as h5file: 
            try: 
                table = h5file.root.Header
            except 'tables.exceptions.NoSuchNodeError':
                print 'This file is in the wrong format'
                return
            
            for h in table.iterrows():
                self.header.unitnum = h['unitnum'][0]
                self.header.version = h['version'][0]
                
                timestamp = h['timestamp'][0]
                
                self.header.datetime = datetime.datetime.fromtimestamp(timestamp)            
        
                self.header.shotsum = h['shotsum'][0]  #total number of shots collected per profile
                self.header.trigfreq = h['trigfreq'][0] #laser trigger frequency (usually 2500 Hz)
                self.header.energy = h['energy'][0] #mean of laser energy monitor in uJ
                self.header.temp_0 = h['temp_0'][0]  #mean of A/D#0 readings*100
                self.header.temp_1 = h['temp_1'][0]  #mean of A/D#1 readings*100
                self.header.temp_2 = h['temp_2'][0]  #mean of A/D#2 readings*100
                self.header.temp_3 = h['temp_3'][0]  #mean of A/D#3 readings*100
                self.header.temp_4 = h['temp_4'][0]  #mean of A/D#4 readings*100
                
                self.header.bg_avg1 = h['bg_avg1'][0] #mean background signal value for channel 1
                self.header.bg_std1 = h['bg_std1'][0] #standard deviation of backgruond signal for channel 1
        
                self.header.numchans = h['numchans'][0] #number of channels
                self.header.numbins = h['numbins'][0] #total number of bins per channel
        
                self.header.bintime = h['bintime'][0]  #bin width in seconds
                
                self.header.rangecal = h['rangecal'][0] #range offset in meters, default is 0
        
                self.header.databins = h['databins'][0]  #number of bins not including those used for background
                self.header.scanflag = h['scanflag'][0]  #0: no scanner, 1: scanner
                self.header.backbins = h['backbins'][0]  #number of background bins
        
                self.header.az = h['az'][0]  #scanner azimuth angle
                self.header.el = h['el'][0]  #scanner elevation angle
                self.header.deg = h['deg'][0] #compass degrees (currently unused)
                self.header.pvolt0 = h['pvolt0'][0] #currently unused
                self.header.pvolt1 = h['pvolt1'][0] #currently unused
                self.header.gpslat = h['gpslat'][0] #GPS latitude in decimal degreees (-999.0 if no GPS)
                self.header.gpslon = h['gpslon'][0] #GPS longitude in decimal degrees (-999.0 if no GPS)
                self.header.cloudbase = h['cloudbase'][0] #cloud base height in [m]
        
                self.header.baddat = h['baddat'][0]  #0: good data, 1: bad data
                self.header.version = h['version'][0] #version of file format.  current version is 1
        
                self.header.bg_avg2 = h['bg_avg2'][0] #mean background signal for channel 2
                self.header.bg_std2 = h['bg_std2'][0] #mean background standard deviation for channel 2
        
                self.header.mcs = h['mcs'][0]  #MCS mode register  Bit#7: 0-normal, 1-polarization
                                             #Bit#6-5: polarization toggling: 00-linear polarizer control
                                             #01-toggling pol control, 10-toggling pol control 11-circular pol control
        
                self.header.firstbin = h['firstbin'][0]  #bin # of first return data
                self.header.systype = h['systype'][0]  #0: standard MPL, 1: mini MPL
                self.header.syncrate = h['syncrate'][0]  #mini-MPL only, sync pulses seen per second
                self.header.firstback = h['firstback'][0] #mini-MPL only, first bin used for background calcs
                self.header.headersize2 = h['headersize2'][0] #size of additional header data (currently unused)
     
        copoldat = pan.read_hdf(filename,'copol')
        crosspoldat = pan.read_hdf(filename,'crosspol')
            
        self.data = [copoldat,crosspoldat]
        
        return self
    
    def fromMPL(self, filename):
        
        binfile = open(filename,'rb')
        
        intarray16 = array.array('H')
        intarray32 = array.array('L')
        floatarray = array.array('f')
        byte_array = array.array('B')
        datavals = array.array('f')
        
        
        intarray16.fromfile(binfile, 8)
        intarray32.fromfile(binfile, 8)
        floatarray.fromfile(binfile, 2)
        intarray16.fromfile(binfile, 1)
        intarray32.fromfile(binfile, 1)
        floatarray.fromfile(binfile, 2)
        intarray16.fromfile(binfile, 3)
        floatarray.fromfile(binfile, 8)
        byte_array.fromfile(binfile, 2)
        floatarray.fromfile(binfile, 2)
        byte_array.fromfile(binfile, 1)
        intarray16.fromfile(binfile, 1)
        byte_array.fromfile(binfile, 1)
        intarray16.fromfile(binfile, 3)
                   
        self.header.unitnum = intarray16[0]
        self.header.version = intarray16[1]
        year = intarray16[2]
        month = intarray16[3]
        day = intarray16[4]
        hour = intarray16[5]
        minute = intarray16[6]
        second = intarray16[7]

        self.header.datetime = datetime.datetime(year,month,day,hour,minute,second)

        self.header.shotsum = intarray32[0]  #total number of shots collected per profile
        self.header.trigfreq = intarray32[1] #laser trigger frequency (usually 2500 Hz)
        self.header.energy = intarray32[2]/1000.0  #mean of laser energy monitor in uJ                      
        self.header.temp_0 = intarray32[3]/1000.0  #mean of A/D#0 readings*100
        self.header.temp_1 = intarray32[4]/1000.0  #mean of A/D#1 readings*100
        self.header.temp_2 = intarray32[5]/1000.0  #mean of A/D#2 readings*100
        self.header.temp_3 = intarray32[6]/1000.0  #mean of A/D#3 readings*100
        self.header.temp_4 = intarray32[7]/1000.0  #mean of A/D#4 readings*100
        
        self.header.bg_avg1 = floatarray[0] #mean background signal value for channel 1
        self.header.bg_std1 = floatarray[1] #standard deviation of backgruond signal for channel 1

        self.header.numchans = intarray16[8] #number of channels
        self.header.numbins = intarray32[8] #total number of bins per channel

        self.header.bintime = floatarray[2]  #bin width in seconds
        
        self.header.rangecal = floatarray[3] #range offset in meters, default is 0

        self.header.databins = intarray16[9]  #number of bins not including those used for background
        self.header.scanflag = intarray16[10]  #0: no scanner, 1: scanner
        self.header.backbins = intarray16[11]  #number of background bins

        self.header.az = floatarray[4]  #scanner azimuth angle
        self.header.el = floatarray[5]  #scanner elevation angle
        self.header.deg = floatarray[6] #compass degrees (currently unused)
        self.header.pvolt0 = floatarray[7] #currently unused
        self.header.pvolt1 = floatarray[8] #currently unused
        self.header.gpslat = floatarray[9] #GPS latitude in decimal degreees (-999.0 if no GPS)
        self.header.gpslon = floatarray[10]#GPS longitude in decimal degrees (-999.0 if no GPS)
        self.header.cloudbase = floatarray[11] #cloud base height in [m]

        self.header.baddat = byte_array[0]  #0: good data, 1: bad data
        self.header.version = byte_array[1] #version of file format.  current version is 1

        self.header.bg_avg2 = floatarray[12] #mean background signal for channel 2
        self.header.bg_std2 = floatarray[13] #mean background standard deviation for channel 2

        self.header.mcs = byte_array[2]  #MCS mode register  Bit#7: 0-normal, 1-polarization
                                     #Bit#6-5: polarization toggling: 00-linear polarizer control
                                     #01-toggling pol control, 10-toggling pol control 11-circular pol control

        self.header.firstbin = intarray16[12]  #bin # of first return data
        self.header.systype = byte_array[3]   #0: standard MPL, 1: mini MPL
        self.header.syncrate = intarray16[13]  #mini-MPL only, sync pulses seen per second
        self.header.firstback = intarray16[14] #mini-MPL only, first bin used for background calcs
        self.header.headersize2 = intarray16[15] #size of additional header data (currently unused)
                 
        numbins = self.header.numbins
        numchans = self.header.numchans          
        
        while True:
            try:
                datavals.fromfile(binfile,(numbins*numchans))                    
            except EOFError:
                break
            
        altstep = self.header.bintime*const.c #altitude step in meters
        maxalt = numbins*altstep
        minalt = self.header.rangecal
        altrange = np.arange(minalt,maxalt,altstep,dtype='float')
        
        timestep = datetime.timedelta(seconds = self.header.shotsum/self.header.trigfreq) #time between profiles in seconds
        dt = self.header.datetime
      
        if numchans == 2:
            
            profdat_copol = []
            profdat_crosspol = []
            
            n = 0
            while n < len(datavals):
                profdat_crosspol.append(datavals[n:n+numbins].tolist()) 
                n += numbins                               
                profdat_copol.append(datavals[n:n+numbins].tolist())
                
                n += numbins + 32
                
                #note: this more adaptive algorithm should have worked
                #I'm not sure why it doesn't so I'm going with the
                #patch for now - 04-25-13
#                            m = n+numbins
#                            temp_prof_crosspol = datarray[n:m]
#                            temp_prof_copol = datarray[m:m+numbins]
#                            
#                            
#                            if all(temp_prof_crosspol > 1e-6):                             
#                                profdat_crosspol.append(temp_prof_crosspol)                                
#                                profdat_copol.append(temp_prof_copol)
#                                n += 2*numbins
#                            else:
#                                n += 1
            copoldat = np.array(profdat_copol)
            crosspoldat = np.array(profdat_crosspol)
            numsteps = np.shape(copoldat)[0]

            timerange = []
            for t in range(0,numsteps):
                timerange.append(dt)
                dt += timestep  
                
            self.data = [pan.DataFrame(copoldat,index = timerange,columns = altrange),\
            pan.DataFrame(crosspoldat,index = timerange,columns = altrange)]
            
        else:
            raise ValueError('Wrong number of channels')
        
        return self

    
    def copy(self):
        
        MPLout = MPL()
        MPLout = deepcopy(self)
        
        return MPLout
    
    def append(self,MPLnew):
        
        if not self.header:
            self.header = MPLnew.header
        for n in range(self.header.numchans):
            self.data[n] = self.data[n].append(MPLnew.data[n])
        
        return self
    
    def save_to_MPL(self,filename):
        
        with open(filename, 'wb') as MPLout:
        
            intarray16 = array.array('H')
            intarray32 = array.array('L')
            floatarray = array.array('f')
            byte_array = array.array('B')
            datavals = array.array('f')
                       
            intarray16.append(self.header.unitnum)
            intarray16.append(self.header.version)
            d = self.header.datetime
            
            intarray16.append(d.year)
            intarray16.append(d.month)
            intarray16.append(d.day)
            intarray16.append(d.hour)
            intarray16.append(d.minute)
            intarray16.append(d.second)
    
            intarray32.append(self.header.shotsum) #total number of shots collected per profile
            intarray32.append(self.header.trigfreq) #laser trigger frequency (usually 2500 Hz)
            intarray32.append(int(self.header.energy*1000))  #mean of laser energy monitor in uJ                      
            intarray32.append(int(self.header.energy*1000))  #mean of A/D#0 readings*100
            intarray32.append(int(self.header.temp_1*1000)) #mean of A/D#1 readings*100
            intarray32.append(int(self.header.temp_2*1000))  #mean of A/D#2 readings*100
            intarray32.append(int(self.header.temp_3*1000))  #mean of A/D#3 readings*100
            intarray32.append(int(self.header.temp_4*1000))  #mean of A/D#4 readings*100
            
            floatarray.append(self.header.bg_avg1) #mean background signal value for channel 1
            floatarray.append(self.header.bg_std1) #standard deviation of backgruond signal for channel 1
    
            intarray16.append(self.header.numchans) #number of channels
            intarray32.append(self.header.numbins) #total number of bins per channel
    
            floatarray.append(self.header.bintime)  #bin width in seconds
            
            floatarray.append(self.header.rangecal) #range offset in meters, default is 0
    
            intarray16.append(self.header.databins) #number of bins not including those used for background
            intarray16.append(self.header.scanflag)  #0: no scanner, 1: scanner
            intarray16.append(self.header.backbins) #number of background bins
    
            floatarray.append(self.header.az)  #scanner azimuth angle
            floatarray.append(self.header.el) #scanner elevation angle
            floatarray.append(self.header.deg) #compass degrees (currently unused)
            floatarray.append(self.header.pvolt0) #currently unused
            floatarray.append(self.header.pvolt1) #currently unused
            floatarray.append(self.header.gpslat) #GPS latitude in decimal degreees (-999.0 if no GPS)
            floatarray.append(self.header.gpslon) #GPS longitude in decimal degrees (-999.0 if no GPS)
            floatarray.append(self.header.cloudbase) #cloud base height in [m]
    
            byte_array.append(self.header.baddat)  #0: good data, 1: bad data
            byte_array.append(self.header.version) #version of file format.  current version is 1
    
            floatarray.append(self.header.bg_avg2) #mean background signal for channel 2
            floatarray.append(self.header.bg_std2) #mean background standard deviation for channel 2
    
            byte_array.append(self.header.mcs) #MCS mode register  Bit#7: 0-normal, 1-polarization
                                         #Bit#6-5: polarization toggling: 00-linear polarizer control
                                         #01-toggling pol control, 10-toggling pol control 11-circular pol control
    
            intarray16.append(self.header.firstbin)  #bin # of first return data
            byte_array.append(self.header.systype)  #0: standard MPL, 1: mini MPL
            intarray16.append(self.header.syncrate)  #mini-MPL only, sync pulses seen per second
            intarray16.append(self.header.firstback) #mini-MPL only, first bin used for background calcs
            intarray16.append(self.header.headersize2) #size of additional header data (currently unused)
                     
            copoldat = self.data[0].values
            crosspoldat = self.data[1].values
            
            profile_buffer = np.zeros(32,dtype='float')
            
            for n in range(np.shape(copoldat)[0]):
                datavals.fromlist(crosspoldat[n].tolist())
                datavals.fromlist(copoldat[n].tolist())
                datavals.fromlist(profile_buffer.tolist())
        
            temparray = intarray16[:8]                                    
            temparray.tofile(MPLout)
            temparray = intarray32[:8]                 
            temparray.tofile(MPLout)
            temparray = floatarray[:2]             
            temparray.tofile(MPLout)
            temparray = array.array('H',[intarray16[8]])             
            temparray.tofile(MPLout)
            temparray = array.array('L',[intarray32[8]])             
            temparray.tofile(MPLout)
            temparray = floatarray[2:4]            
            temparray.tofile(MPLout)
            temparray = intarray16[9:12]             
            temparray.tofile(MPLout)
            temparray = floatarray[4:12]
            temparray.tofile(MPLout)
            temparray = byte_array[:2]            
            temparray.tofile(MPLout)
            temparray = floatarray[12:14]            
            temparray.tofile(MPLout)        
            temparray = array.array('B',[byte_array[2]])             
            temparray.tofile(MPLout)
            temparray = array.array('H',[intarray16[12]])             
            temparray.tofile(MPLout)
            temparray = array.array('B',[byte_array[3]])             
            temparray.tofile(MPLout)
            temparray = intarray16[13:]             
            temparray.tofile(MPLout)
            datavals.tofile(MPLout)
        
    def save_to_HDF(self, filename, appendflag = 'false'):
        
        #Define class header to create columns to hold header data from .mpl binary file
        class header(tables.IsDescription):
            unitnum = tables.UInt32Col(1)
            version = tables.UInt32Col(1)
            
            timestamp = tables.Time32Col(1)
            
            shotsum = tables.UInt32Col(1)  #total number of shots collected per profile
            trigfreq = tables.UInt32Col(1) #laser trigger frequency (usually 2500 Hz)
            energy = tables.UInt32Col(1) #mean of laser energy monitor in uJ
            temp_0 = tables.UInt32Col(1)  #mean of A/D#0 readings*100
            temp_1 = tables.UInt32Col(1)  #mean of A/D#1 readings*100
            temp_2 = tables.UInt32Col(1)  #mean of A/D#2 readings*100
            temp_3 = tables.UInt32Col(1)  #mean of A/D#3 readings*100
            temp_4 = tables.UInt32Col(1)  #mean of A/D#4 readings*100
            
            bg_avg1 = tables.Float32Col(1) #mean background signal value for channel 1
            bg_std1 = tables.Float32Col(1) #standard deviation of backgruond signal for channel 1
    
            numchans = tables.UInt16Col(1) #number of channels
            numbins = tables.UInt32Col(1) #total number of bins per channel
    
            bintime = tables.Float32Col(1)  #bin width in seconds
            
            rangecal = tables.Float32Col(1) #range offset in meters, default is 0
    
            databins = tables.UInt16Col(1)  #number of bins not including those used for background
            scanflag = tables.UInt16Col(1)  #0: no scanner, 1: scanner
            backbins = tables.UInt16Col(1)  #number of background bins
    
            az = tables.Float32Col(1)  #scanner azimuth angle
            el = tables.Float32Col(1)  #scanner elevation angle
            deg = tables.Float32Col(1) #compass degrees (currently unused)
            pvolt0 = tables.Float32Col(1) #currently unused
            pvolt1 = tables.Float32Col(1) #currently unused
            gpslat = tables.Float32Col(1) #GPS latitude in decimal degreees (-999.0 if no GPS)
            gpslon = tables.Float32Col(1) #GPS longitude in decimal degrees (-999.0 if no GPS)
            cloudbase = tables.Float32Col(1) #cloud base height in [m]
    
            baddat = tables.BoolCol(1)  #0: good data, 1: bad data
            version = tables.BoolCol(1) #version of file format.  current version is 1
    
            bg_avg2 = tables.Float32Col(1) #mean background signal for channel 2
            bg_std2 = tables.Float32Col(1) #mean background standard deviation for channel 2
    
            mcs = tables.BoolCol(8)  #MCS mode register  Bit#7: 0-normal, 1-polarization
                                         #Bit#6-5: polarization toggling: 00-linear polarizer control
                                         #01-toggling pol control, 10-toggling pol control 11-circular pol control
    
            firstbin = tables.UInt16Col(1)  #bin # of first return data
            systype = tables.BoolCol(1)   #0: standard MPL, 1: mini MPL
            syncrate = tables.UInt16Col(1)  #mini-MPL only, sync pulses seen per second
            firstback = tables.UInt16Col(1) #mini-MPL only, first bin used for background calcs
            headersize2 = tables.UInt16Col(1) #size of additional header data (currently unused)
        
        
        with tables.openFile(filename, mode = 'w', title = 'MPL data file') as h5filename:
            
            headertbl = h5filename.create_table('/','Header',header,'Ancillary Data')
              
            headerdat = headertbl.row
                      
            headerdat['unitnum'] = self.header.unitnum
            headerdat['version'] = self.header.version       
            headerdat['timestamp'] = time.mktime(self.header.datetime.timetuple())
        
            headerdat['shotsum'] = self.header.shotsum  #total number of shots collected per profile
            headerdat['trigfreq'] = self.header.trigfreq #laser trigger frequency (usually 2500 Hz)
            headerdat['energy'] = self.header.energy  #mean of laser energy monitor in uJ
            headerdat['temp_0'] = self.header.temp_0  #mean of A/D#0 readings*100
            headerdat['temp_1'] = self.header.temp_1  #mean of A/D#1 readings*100
            headerdat['temp_2'] = self.header.temp_2  #mean of A/D#2 readings*100
            headerdat['temp_3'] = self.header.temp_3  #mean of A/D#3 readings*100
            headerdat['temp_4'] = self.header.temp_4  #mean of A/D#4 readings*100
            
            headerdat['bg_avg1'] = self.header.bg_avg1 #mean background signal value for channel 1
            headerdat['bg_std1'] = self.header.bg_std1 #standard deviation of backgruond signal for channel 1
        
            headerdat['numchans'] = self.header.numchans #number of channels
            headerdat['numbins'] = self.header.numbins #total number of bins per channel
        
            headerdat['bintime'] = self.header.bintime  #bin width in seconds
            
            headerdat['rangecal'] = self.header.rangecal #range offset in meters, default is 0
        
            headerdat['databins'] = self.header.databins  #number of bins not including those used for background
            headerdat['scanflag'] = self.header.scanflag  #0: no scanner, 1: scanner
            headerdat['backbins'] = self.header.backbins  #number of background bins
        
            headerdat['az'] = self.header.az  #scanner azimuth angle
            headerdat['el'] = self.header.el  #scanner elevation angle
            headerdat['deg'] = self.header.deg #compass degrees (currently unused)
            headerdat['pvolt0'] = self.header.pvolt0 #currently unused
            headerdat['pvolt1'] = self.header.pvolt1 #currently unused
            headerdat['gpslat'] = self.header.gpslat #GPS latitude in decimal degreees (-999.0 if no GPS)
            headerdat['gpslon'] = self.header.gpslon #GPS longitude in decimal degrees (-999.0 if no GPS)
            headerdat['cloudbase'] = self.header.cloudbase #cloud base height in [m]
        
            headerdat['baddat'] = self.header.baddat  #0: good data, 1: bad data
            headerdat['version'] = self.header.version #version of file format.  current version is 1
        
            headerdat['bg_avg2'] = self.header.bg_avg2 #mean background signal for channel 2
            headerdat['bg_std2'] = self.header.bg_std2 #mean background standard deviation for channel 2
        
            headerdat['mcs'] = self.header.mcs  #MCS mode register  Bit#7: 0-normal, 1-polarization
                                         #Bit#6-5: polarization toggling: 00-linear polarizer control
                                         #01-toggling pol control, 10-toggling pol control 11-circular pol control
        
            headerdat['firstbin'] = self.header.firstbin  #bin # of first return data
            headerdat['systype'] = self.header.systype   #0: standard MPL, 1: mini MPL
            headerdat['syncrate'] = self.header.syncrate  #mini-MPL only, sync pulses seen per second
            headerdat['firstback'] = self.header.firstback #mini-MPL only, first bin used for background calcs
            headerdat['headersize2'] = self.header.headersize2 #size of additional header data (currently unused)
            
            headerdat.append()
            headertbl.flush()
            
        df_copol = self.data[0]
        df_crosspol = self.data[1]
            
        store = pan.HDFStore(filename)
        store['copol'] = df_copol
        store['crosspol'] = df_crosspol  
        store.close()

    
    def quicklook(self, vlim = []):    
        
        copol = self.data[0]
        crosspol = self.data[1]
        
        print 'Copol Shape is'
        print np.shape(copol)
        print 'Crosspol Shape is'
        print np.shape(crosspol)
        
        cmap = cm.jet
        
        if vlim:
            cmap.set_over('r')
            cmap.set_under('k')
            vmin= vlim[0]
            vmax= vlim[1]
            the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
            fig = plt.figure()
            ax1 = fig.add_subplot(121)
            im1 = ax1.pcolormesh(copol.T,cmap=cmap,norm=the_norm)
            cb = fig.colorbar(im1,extend='both')
        
            ax2 = fig.add_subplot(122)
            im2 = ax2.pcolormesh(crosspol.T,cmap=cmap,norm=the_norm)
            cb = fig.colorbar(im2,extend='both')
        else:
            fig = plt.figure()
            ax1 = fig.add_subplot(121)
            im1 = ax1.pcolormesh(copol.T,cmap=cmap)
            cb = fig.colorbar(im1,extend='both')
        
            ax2 = fig.add_subplot(122)
            im2 = ax2.pcolormesh(crosspol.T,cmap=cmap)
            cb = fig.colorbar(im2,extend='both')
        
        plt.show()
    
    def alt_resample(self, altrange):
        #takes a pandas dataframe generated by mplreader and resamples on regular
        #intervals in altitude and resets the limits of the set
        #note: limits of altrange must be within original limits of altitude data
        
        dataout = []
        for n in range(self.header.numchans):
            print 'Altitude step resampling in progress ...'
            df = self.data[n]
            x = df.columns
            
            numrows = np.size(df.index)
            numcols = np.size(altrange)
        
            newvalues = np.empty([numrows, numcols])
            r = 0
        
            for row in df.iterrows():
                f = interp1d(x,row[1].values)
                newvalues[r,:] = f(altrange)
                r += 1
            dataout.append(pan.DataFrame(data = newvalues, index = df.index, columns = altrange))
            print '... Done!'
            
        self.data = dataout
        self.header.databins = len(altrange)
        self.header.numbins = self.header.databins + self.header.backbins
        self.header.bintime = (altrange[1]-altrange[0])/const.c
        
        return self
    
    def time_resample(self, timestep, timerange = False, s_mode = 'mean'):
        #resamples a pandas dataframe generated by mplreader on a regular timestep
        #and optionally limits it to a preset time range
        #timestep must be in timeseries period format: numF where num=step size and
        #F = offset alias.  Ex: H = hours, M = minutes, S = seconds, L = millieconds
        
        dataout = np.empty(2, dtype=object)
        for n in range(self.header.numchans):    
            print 'Time step regularization in progress ...'
            df = self.data[n]
            if timerange:
                start_time = timerange[0]
                end_time = timerange[1]
        
                dfout = df[(df.index>=start_time) & (df.index<=end_time)]
        
            dataout[n] = dfout.resample(timestep, how = 'mean')
        
            print '... Done!'
            
        self.data = dataout
        df = dataout[0]
        ts = (df.index[2]-df.index[1]).total_seconds()
        
        self.header.shotsum = int(ts*self.header.trigfreq)
        
        return self
        
    def background_subtract(self):
        
        backbins = self.header.backbins
              
        dataout = []
        for n in range(self.header.numchans):
            df = self.data[n]
    
            newvalues = np.empty_like(df.values)
            
            r = 0
            for row in df.iterrows():
                bg = np.mean(row[1].values[-backbins:])
                newvalues[r,:] = np.array(row[1].values - bg)
                r += 1
            
            dataout.append(pan.DataFrame(data = newvalues, index = df.index, columns = df.columns))
        
                
        self.data = dataout
        
        return self
    
    def range_cor(self):
           
        dataout = []
        for n in [0,1]:
            df = self.data[n]
            rsq = np.array(df.columns, dtype=float)**2
            newvalues = np.empty_like(df.values)
            r = 0
            for row in df.iterrows():
                newvalues[r,:] = row[1].values*rsq
                r += 1
            
            dataout.append(pan.DataFrame(data = newvalues, index = df.index, columns = df.columns))
        
        self.data = dataout
        
        return self

    def power_cor(self):
        e = self.header.energy
        f = self.header.trigfreq
        s = self.header.shotsum
        
        
        p = e*f*s  #power in uW = energy per shot [uJ] * shot frequency [Hz] * number of shots per profile
        
        dataout = []
        for n in [0,1]:
            df = self.data[n]
            newvalues = np.empty_like(df.values)
            r = 0
            for row in df.iterrows():
                newvalues[r,:] = row[1].values/p
                r += 1
            
            dataout.append(pan.DataFrame(data = newvalues, index = df.index, columns = df.columns))
        
        self.data = dataout
        
        return self
    
   


if __name__ == '__main__':
    
    import os
    import MPLtools as mtools
    
    olddir = os.getcwd()
    
    os.chdir('C:\Program Files (x86)\SigmaMPL\DATA\\2013_04_25')
    mdat = mtools.MPL()
    mdat.fromMPL('test.mpl')
    
    mdat.save_to_MPL('test.mpl')
    print 'done'
    
    os.chdir(olddir)