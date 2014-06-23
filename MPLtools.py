 # -*- coding: utf-8 -*-
"""
mpltools.py
A bag of tools to be used in processing and interpreting MPL class data
collected by miniMPL
Created on Wed Apr 24 12:08:57 2013

@author: Paul Cottle

"""



def set_dir(titlestring):
    from Tkinter import Tk
    import tkFileDialog
     
    # Make a top-level instance and hide since it is ugly and big.
    root = Tk()
    root.withdraw()
    
    # Make it almost invisible - no decorations, 0 size, top left corner.
    root.overrideredirect(True)
    root.geometry('0x0+0+0')
#    
    # Show window again and lift it to top so it can get focus,
    # otherwise dialogs will end up behind the terminal.
    root.deiconify()
    root.attributes("-topmost",1)
    root.focus_force()
    
    file_path = tkFileDialog.askdirectory(parent=root,title=titlestring)
     
    if file_path != "":
       return str(file_path)
     
    else:
       print "you didn't open anything!"
    
    # Get rid of the top-level instance once to make it actually invisible.
    root.destroy() 
     

       
    
     
def get_files(titlestring,filetype = ('.txt','*.txt')):
    from Tkinter import Tk
    import tkFileDialog
    import re
     
     
    # Make a top-level instance and hide since it is ugly and big.
    root = Tk()
    root.withdraw()
    
    # Make it almost invisible - no decorations, 0 size, top left corner.
    root.overrideredirect(True)
    root.geometry('0x0+0+0')
#    
    # Show window again and lift it to top so it can get focus,
    # otherwise dialogs will end up behind the terminal.
    root.deiconify()
    root.attributes("-topmost",1)
    root.focus_force()

    filenames = []
     
    filenames = tkFileDialog.askopenfilename(title=titlestring, filetypes=[filetype],multiple='True')
    
    #do nothing if already a python list
    if filenames == "": 
        print "You didn't open anything!"  
        return
        
    if isinstance(filenames,list): return filenames

    #http://docs.python.org/library/re.html
    #the re should match: {text and white space in brackets} AND anynonwhitespacetokens
    #*? is a non-greedy match for any character sequence
    #\S is non white space

    #split filenames string up into a proper python list
    result = re.findall("{.*?}|\S+",filenames)

    #remove any {} characters from the start and end of the file names
    result = [ re.sub("^{|}$","",i) for i in result ]     
    return result

    
    root.destroy()
      
    
def MPLtoHDF(filename, appendflag = 'False'):
    import numpy as np
#    import tables
    import array
    import pandas as pan
    import datetime
    from scipy import constants as const
    
    h5filename = filename.split('.')[0]+'_proc.h5'

    with open(filename,'rb') as binfile:
    
        profdat_copol = {}
        profdat_crosspol = {}
        header = {}
        headerdat = {}
        
        profnum = 0
        
        while True:
            try:
                intarray16 = array.array('H')
                intarray32 = array.array('I') # L is 8 byte on Xenon
                floatarray = array.array('f')
                byte_array = array.array('B')
                copolvals = array.array('f')  
                crosspolvals = array.array('f')  
                
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
                
                if headerdat['headersize2'] > 128:
                    byte_array.fromfile(binfile, 1)
                    floatarray.fromfile(binfile, 6)
                    intarray16.fromfile(binfile, 1)
                    floatarray.fromfile(binfile, 2)
                
                    headerdat['Weatherstat'] = byte_array[4]   #0:weatehr station not used, 1:weather station used
                    headerdat['Int_temp'] = floatarray[14]  #Temperature inside in deg. Celsius
                    headerdat['Ext_temp'] = floatarray[15]  #Temp. outside
                    headerdat['Int_humid'] = floatarray[16] #Humidity inside in %
                    headerdat['Ext_humid'] = floatarray[17] #Hum. outside
                    headerdat['Dewpoint'] = floatarray[18]  #dewpoint in deg. Celsius
                    headerdat['Wnd_spd'] = floatarray[19]  #wind speed in km/h
                    headerdat['Wnd_dir'] = intarray16[16]  #wind direction in deg.
                    headerdat['Press'] = floatarray[20]  #Barometric pressure in hPa
                    headerdat['Rain'] = floatarray[21]  #Rain rate in mm/hr
                
                
                numbins = headerdat['numbins']
                numchans = headerdat['numchans'] 
                altstep = headerdat['bintime']*const.c #altitude step in meters
                maxalt = numbins*altstep
                minalt = headerdat['rangecal']
                altrange = np.arange(minalt,maxalt,altstep,dtype='float')
                
                if numchans == 2:
                    copolvals.fromfile(binfile, numbins) 
                    temp = np.array(copolvals)
                    profdat_crosspol[dt] = pan.Series(temp, index = altrange)
                    crosspolvals.fromfile(binfile, numbins) 
                    temp = np.array(crosspolvals)
                    profdat_copol[dt] = pan.Series(temp, index = altrange)
                else:
                    raise ValueError('Wrong number of channels')
                
                headerdat['profnum'] = profnum
                profnum += 1
                header[dt] = pan.Series(headerdat)
            except EOFError:
                break
                
            
        df_copol = pan.DataFrame.from_dict(profdat_copol, orient='index')
        df_crosspol = pan.DataFrame.from_dict(profdat_crosspol, orient='index')
        df_header = pan.DataFrame.from_dict(header, orient='index')
        
    store = pan.HDFStore(h5filename)
    store['copol_raw'] = df_copol
    store['crosspol_raw'] = df_crosspol
    store['header'] = df_header
    store.close()
       

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
        
    def __init__(self,filename=[]):
        
        self.data = [] #slot for lidar raw data array
        self.rsq = []  #slot for range corrected, background subtracted data
        self.NRB = []  #slot for Normalized Relative Backscatter array 
        self.depolrat = []  #slot fo depol ratio array
        self.header = []  #slot for header data
        self.SNR = []
        

    def copy(self):
        #currnetly not working!
        from copy import deepcopy
        
        return MPL(deepcopy(self))
    
    def append(self,MPLnew):
        
        if type(self.header) == 'list':
            self.header = MPLnew.header
        else:
            maxprofnum = self.header['profnum'][-1]
            MPLnew.header['profnum'] = [p+maxprofnum for p in MPLnew.header['profnum']]
            
            self.header = self.header.append(MPLnew.header)
        
        if not self.data:
            self.data = MPLnew.data
        else:            
            for n in range(self.header['numchans'][0]):
                self.data[n] = self.data[n].append(MPLnew.data[n])
        
        if MPLnew.rsq:
            if not self.rsq:
                self.rsq = MPLnew.rsq
            else:
                for n in range(self.header['numchans'][0]):
                    self.rsq[n] = self.rsq[n].append(MPLnew.rsq[n])

        if MPLnew.NRB:
            if not self.NRB:
                self.NRB = MPLnew.NRB
            else:
                for n in range(self.header['numchans'][0]):
                    self.NRB[n] = self.NRB[n].append(MPLnew.NRB[n])
        
        if MPLnew.depolrat:
            if not self.depolrat:
                self.depolrat = MPLnew.depolrat
            else:
                self.depolrat[0] = self.depolrat[0].append(MPLnew.depolrat[0])
                    
        return self
    
    def fromMPL(self, filename):
        import numpy as np
        import datetime
        import array
        import pandas as pan
        from scipy import constants as const
        
        with open(filename,'rb') as binfile:
            profdat_copol = {}
            profdat_crosspol = {}
            header = {}
            
            profnum = 0
            
            while True:
                try:
                    intarray16 = array.array('H')
                    intarray32 = array.array('I') # L is 8 byte on Xenon
                    floatarray = array.array('f')
                    byte_array = array.array('B')
                    copolvals = array.array('f')  
                    crosspolvals = array.array('f')  
                    
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
                    
                    headerdat = {}
                    headerdat['unitnum'] = intarray16[0]
                    headerdat['version'] = intarray16[1]
                    year = intarray16[2]
                    month = intarray16[3]
                    day = intarray16[4]
                    hour = intarray16[5]
                    minute = intarray16[6]
                    second = intarray16[7]
                                       
                    dt = datetime.datetime(year,month,day,hour,minute,second)
            
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
            
                    # print intarray16.itemsize, intarray32.itemsize # L vs I
            
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
                    headerdat['gpslon'] = floatarray[10] #GPS longitude in decimal degrees (-999.0 if no GPS)
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
                    headerdat['profnum'] = profnum
                    profnum += 1
                    
                    if headerdat['headersize2'] > 128:
                        byte_array.fromfile(binfile, 1)
                        floatarray.fromfile(binfile, 6)
                        intarray16.fromfile(binfile, 1)
                        floatarray.fromfile(binfile, 2)
                    
                        headerdat['Weatherstat'] = byte_array[4]   #0:weatehr station not used, 1:weather station used
                        headerdat['Int_temp'] = floatarray[14]  #Temperature inside in deg. Celsius
                        headerdat['Ext_temp'] = floatarray[15]  #Temp. outside
                        headerdat['Int_humid'] = floatarray[16] #Humidity inside in %
                        headerdat['Ext_humid'] = floatarray[17] #Hum. outside
                        headerdat['Dewpoint'] = floatarray[18]  #dewpoint in deg. Celsius
                        headerdat['Wnd_spd'] = floatarray[19]  #wind speed in km/h
                        headerdat['Wnd_dir'] = intarray16[16]  #wind direction in deg.
                        headerdat['Press'] = floatarray[20]  #Barometric pressure in hPa
                        headerdat['Rain'] = floatarray[21]  #Rain rate in mm/hr

                    numbins = headerdat['numbins']
                    numchans = headerdat['numchans'] 
                    altstep = headerdat['bintime']*const.c/2 #altitude step in meters
                    maxalt = numbins*altstep
                    minalt = headerdat['rangecal']
                    altrange = np.arange(minalt,maxalt,altstep,dtype='float')
                    
                    if numchans == 2:
                        crosspolvals.fromfile(binfile, numbins) 
                        temp = np.array(crosspolvals)                       
                        profdat_crosspol[dt] = temp
                        copolvals.fromfile(binfile, numbins) 
                        temp = np.array(copolvals)
                        profdat_copol[dt] = temp
                    else:
                        raise ValueError('Wrong number of channels')
                    
                    header[dt] = headerdat
                
                except EOFError:
                    break
            
        df_copol = pan.DataFrame.from_dict(profdat_copol,orient = 'index')
        df_copol.columns = altrange
        df_crosspol = pan.DataFrame.from_dict(profdat_crosspol,orient = 'index')
        df_crosspol.columns = altrange
        df_header = pan.DataFrame.from_dict(header, orient = 'index')
                
        self.data = [df_copol, df_crosspol]
        
        
        self.header = df_header

        return self        
    
    def fromHDF(self, filename, verbose = False):
        import pandas as pan
                
        copoldat = pan.read_hdf(filename,'copol_raw')
        crosspoldat = pan.read_hdf(filename,'crosspol_raw')
        header = pan.read_hdf(filename,'header')
            
        self.data = [copoldat,crosspoldat]
        self.header = header
        
        try:
            copoldat_rsq = pan.read_hdf(filename,'copol_rsq')
            crosspoldat_rsq = pan.read_hdf(filename,'crosspol_rsq')
            self.rsq = [copoldat_rsq,crosspoldat_rsq]
        except KeyError:
            if verbose:
                print "Warning:  No Range-squared file"
        
        try:
            copoldat_NRB = pan.read_hdf(filename,'copol_NRB')
            crosspoldat_NRB = pan.read_hdf(filename,'crosspol_NRB')
            self.NRB = [copoldat_NRB,crosspoldat_NRB]
        except KeyError:
            if verbose:
                print "Warning: No NRB file"
        
        try:
            self.depolrat = [pan.read_hdf(filename,'depolrat')]
        except KeyError:
            if verbose:
                print "Warning: No Depol Ratio file"

        return self
    
#    def save_to_MPL(self,filename):
#        #currently not working!
#        import numpy as np
#        import array
#        
#        with open(filename, 'wb') as MPLout:
#        
#            intarray16 = array.array('H')
#            intarray32 = array.array('L')
#            floatarray = array.array('f')
#            byte_array = array.array('B')
#            datavals = array.array('f')
#                       
#            intarray16.append(self.header.unitnum)
#            intarray16.append(self.header.version)
#            d = self.header.datetime
#            
#            intarray16.append(d.year)
#            intarray16.append(d.month)
#            intarray16.append(d.day)
#            intarray16.append(d.hour)
#            intarray16.append(d.minute)
#            intarray16.append(d.second)
#    
#            intarray32.append(self.header.shotsum) #total number of shots collected per profile
#            intarray32.append(self.header.trigfreq) #laser trigger frequency (usually 2500 Hz)
#            intarray32.append(int(self.header.energy*1000))  #mean of laser energy monitor in uJ                      
#            intarray32.append(int(self.header.energy*1000))  #mean of A/D#0 readings*100
#            intarray32.append(int(self.header.temp_1*1000)) #mean of A/D#1 readings*100
#            intarray32.append(int(self.header.temp_2*1000))  #mean of A/D#2 readings*100
#            intarray32.append(int(self.header.temp_3*1000))  #mean of A/D#3 readings*100
#            intarray32.append(int(self.header.temp_4*1000))  #mean of A/D#4 readings*100
#            
#            floatarray.append(self.header.bg_avg1) #mean background signal value for channel 1
#            floatarray.append(self.header.bg_std1) #standard deviation of backgruond signal for channel 1
#    
#            intarray16.append(self.header.numchans) #number of channels
#            intarray32.append(self.header.numbins) #total number of bins per channel
#    
#            floatarray.append(self.header.bintime)  #bin width in seconds
#            
#            floatarray.append(self.header.rangecal) #range offset in meters, default is 0
#    
#            intarray16.append(self.header.databins) #number of bins not including those used for background
#            intarray16.append(self.header.scanflag)  #0: no scanner, 1: scanner
#            intarray16.append(self.header.backbins) #number of background bins
#    
#            floatarray.append(self.header.az)  #scanner azimuth angle
#            floatarray.append(self.header.el) #scanner elevation angle
#            floatarray.append(self.header.deg) #compass degrees (currently unused)
#            floatarray.append(self.header.pvolt0) #currently unused
#            floatarray.append(self.header.pvolt1) #currently unused
#            floatarray.append(self.header.gpslat) #GPS latitude in decimal degreees (-999.0 if no GPS)
#            floatarray.append(self.header.gpslon) #GPS longitude in decimal degrees (-999.0 if no GPS)
#            floatarray.append(self.header.cloudbase) #cloud base height in [m]
#    
#            byte_array.append(self.header.baddat)  #0: good data, 1: bad data
#            byte_array.append(self.header.version) #version of file format.  current version is 1
#    
#            floatarray.append(self.header.bg_avg2) #mean background signal for channel 2
#            floatarray.append(self.header.bg_std2) #mean background standard deviation for channel 2
#    
#            byte_array.append(self.header.mcs) #MCS mode register  Bit#7: 0-normal, 1-polarization
#                                         #Bit#6-5: polarization toggling: 00-linear polarizer control
#                                         #01-toggling pol control, 10-toggling pol control 11-circular pol control
#    
#            intarray16.append(self.header.firstbin)  #bin # of first return data
#            byte_array.append(self.header.systype)  #0: standard MPL, 1: mini MPL
#            intarray16.append(self.header.syncrate)  #mini-MPL only, sync pulses seen per second
#            intarray16.append(self.header.firstback) #mini-MPL only, first bin used for background calcs
#            intarray16.append(self.header.headersize2) #size of additional header data (currently unused)
#                     
#            copoldat = self.data[0].values
#            crosspoldat = self.data[1].values
#            
#            profile_buffer = np.zeros(32,dtype='float')
#            
#            for n in range(np.shape(copoldat)[0]):
#                datavals.fromlist(crosspoldat[n].tolist())
#                datavals.fromlist(copoldat[n].tolist())
#                datavals.fromlist(profile_buffer.tolist())
#        
#            temparray = intarray16[:8]                                    
#            temparray.tofile(MPLout)
#            temparray = intarray32[:8]                 
#            temparray.tofile(MPLout)
#            temparray = floatarray[:2]             
#            temparray.tofile(MPLout)
#            temparray = array.array('H',[intarray16[8]])             
#            temparray.tofile(MPLout)
#            temparray = array.array('L',[intarray32[8]])             
#            temparray.tofile(MPLout)
#            temparray = floatarray[2:4]            
#            temparray.tofile(MPLout)
#            temparray = intarray16[9:12]             
#            temparray.tofile(MPLout)
#            temparray = floatarray[4:12]
#            temparray.tofile(MPLout)
#            temparray = byte_array[:2]            
#            temparray.tofile(MPLout)
#            temparray = floatarray[12:14]            
#            temparray.tofile(MPLout)        
#            temparray = array.array('B',[byte_array[2]])             
#            temparray.tofile(MPLout)
#            temparray = array.array('H',[intarray16[12]])             
#            temparray.tofile(MPLout)
#            temparray = array.array('B',[byte_array[3]])             
#            temparray.tofile(MPLout)
#            temparray = intarray16[13:]             
#            temparray.tofile(MPLout)
#            datavals.tofile(MPLout)
        
    def save_to_HDF(self, filename, appendflag = 'false'):
        import pandas as pan
        
        store = pan.HDFStore(filename)
        
        store['header'] = self.header
        
        if self.data:
            df_copol = self.data[0]
            df_crosspol = self.data[1]
            store['copol_raw'] = df_copol
            store['crosspol_raw'] = df_crosspol
        
        if self.rsq:
            df_copol = self.rsq[0]
            df_crosspol = self.rsq[1]
            store['copol_rsq'] = df_copol
            store['crosspol_rsq'] = df_crosspol
        
        if self.NRB:
            df_copol = self.NRB[0]
            df_crosspol = self.NRB[1]
            store['copol_NRB'] = df_copol
            store['crosspol_NRB'] = df_crosspol
        
        if self.depolrat:
            df_depolrat = self.depolrat[0]
            store['depolrat'] = df_depolrat
        
        store.close()

    
    def alt_resample(self,altrange,underfill=True,verbose=False):
        #takes a pandas dataframe generated by mplreader and resamples on regular
        #intervals in altitude and resets the limits of the set
        #note: limits of altrange must be within original limits of altitude data
        import numpy as np
        import pandas as pan
        from scipy import constants as const
        from scipy.interpolate import interp1d
        
        dataout = []
        rsqout = []
        nrbout = []
        depolout=[]
        
        if verbose:
            print 'Altitude step resampling in progress ...'
        for n in range(self.header['numchans'][0]):            
            df = self.data[n]
            x = df.columns
            
            minalt = df.columns[0]
            maxalt = df.columns[-1]
            
            if minalt > altrange[0]:
                altrange = altrange[altrange >= minalt]
                if verbose:
                    print "WARNING: Minimum altitude reset to {0}".format(altrange[0])
            
            if maxalt < altrange[-1]:
                altrange = altrange[altrange <= maxalt]
                if verbose:
                    print "WARNING: Maximum altitude reset to {0}".format(altrange[-1])
                    
            numrows = np.size(df.index)
            numcols = np.size(altrange)
        
            newvalues = np.empty([numrows, numcols])
            r = 0
        
            for row in df.iterrows():
                f = interp1d(x,row[1].values)
                newvalues[r,:] = f(altrange)
                r += 1
            dataout.append(pan.DataFrame(data = newvalues, index = df.index, columns = altrange))
        
        self.data = dataout
        
        if self.rsq:     
            for n in range(self.header['numchans'][0]): 
                df = self.rsq[n]
                x = df.columns
                
                minalt = df.columns[0]
                maxalt = df.columns[-1]
                
                if minalt > altrange[0]:
                    altrange = altrange[altrange >= minalt]
                    if verbose:
                        print "WARNING: Minimum altitude reset to {0}".format(altrange[0])
                
                if maxalt < altrange[-1]:
                    altrange = altrange[altrange <= maxalt]
                    if verbose:
                        print "WARNING: Maximum altitude reset to {0}".format(altrange[-1])
                        
                numrows = np.size(df.index)
                numcols = np.size(altrange)
            
                newvalues = np.empty([numrows, numcols])
                r = 0
            
                for row in df.iterrows():
                    f = interp1d(x,row[1].values)
                    newvalues[r,:] = f(altrange)
                    r += 1
                rsqout.append(pan.DataFrame(data = newvalues, index = df.index, columns = altrange))
            
            self.rsq=rsqout
        else:
            if verbose:
                print "No Range-Squared Profiles"
        
        if self.NRB:       
            for n in range(self.header['numchans'][0]): 
                df = self.NRB[n]
                x = df.columns
                
                minalt = df.columns[0]
                maxalt = df.columns[-1]
                
                if minalt > altrange[0]:
                    altrange = altrange[altrange >= minalt]
                    if verbose:
                        print "WARNING: Minimum altitude reset to {0}".format(altrange[0])
                
                if maxalt < altrange[-1]:
                    altrange = altrange[altrange <= maxalt]
                    if verbose:
                        print "WARNING: Maximum altitude reset to {0}".format(altrange[-1])
                        
                numrows = np.size(df.index)
                numcols = np.size(altrange)
            
                newvalues = np.empty([numrows, numcols])
                r = 0
            
                for row in df.iterrows():
                    f = interp1d(x,row[1].values)
                    newvalues[r,:] = f(altrange)
                    r += 1
                nrbout.append(pan.DataFrame(data = newvalues, index = df.index, columns = altrange))
            
            self.NRB=nrbout
        else:
            if verbose:
                print "No NRB Profiles"    

        if self.depolrat:       
            df = self.depolrat[0]
            x = df.columns
            
            minalt = df.columns[0]
            maxalt = df.columns[-1]
            
            if minalt > altrange[0]:
                altrange = altrange[altrange >= minalt]
                if verbose:
                    print "WARNING: Minimum altitude reset to {0}".format(altrange[0])
            
            if maxalt < altrange[-1]:
                altrange = altrange[altrange <= maxalt]
                if verbose:
                    print "WARNING: Maximum altitude reset to {0}".format(altrange[-1])
                    
            numrows = np.size(df.index)
            numcols = np.size(altrange)
        
            newvalues = np.empty([numrows, numcols])
            r = 0
        
            for row in df.iterrows():
                f = interp1d(x,row[1].values)
                newvalues[r,:] = f(altrange)
                r += 1
            depolout.append(pan.DataFrame(data = newvalues, index = df.index, columns = altrange))
            
            self.depolrat=depolout
        else:
            if verbose:
                print "No Depol Ratio Profiles"  
                
        if verbose:
            print '... Done!'
            
        
        self.header['numbins'] = [len(altrange) for db in self.header['numbins']]
        self.header['databins'] = [db - bb for (db,bb) in zip(self.header['numbins'],self.header['backbins'])]
        self.header['bintime'] = [(altrange[1]-altrange[0])/const.c for bt in self.header['bintime']]
        self.header['firstback'] = [len(altrange)+1 for fb in self.header['firstback']]
        
        return self
    
    def time_resample(self, timestep=[], starttime=[],endtime=[], datamethod = 'mean',verbose=False):
        #resamples a pandas dataframe generated by mplreader on a regular timestep
        #and optionally limits it to a preset time range
        #timestep must be in timeSeries period format: numF where num=step size and
        #F = offset alias.  Ex: H = hours, M = minutes, S = seconds, L = millieconds
        import pandas as pan
        import numpy as np
        
        temphead = {}
        self.header.sort_index(inplace=True)
        if starttime:
            self.header = self.header.loc[self.header.index>=starttime]
        if endtime:
            self.header = self.header.loc[self.header.index<=endtime]
        if timestep:
            for col in self.header:           
                sumcols = ['shotsum']
                firstcols = ['unitnum','version','numchans','scanflag','version','mcs','systype','numchans']
                intmeancols = ['numbins','databins','backbins','firstbin','firstback']
                maxcols = ['baddat']            
                if col in sumcols: headermethod = 'sum'            
                elif col in firstcols: headermethod = 'first'            
                elif col in maxcols:  headermethod = 'max'
                else: headermethod = 'mean'
                            
                temphead[col] = self.header[col].resample(timestep, how = headermethod)            
                
                if (col in intmeancols) or (col in firstcols) or (col in maxcols):
                    try:
                        temphead[col] = temphead[col].astype('int32')
                    except ValueError:
                        whereisna = np.isnan(temphead[col])
                        temphead[col][whereisna] = -999
                        temphead[col] = temphead[col].astype('int32')

        
            self.header = pan.DataFrame(temphead)
        if verbose:
                print 'Time step regularization in progress ...'
                
        for n in range(self.header['numchans'][0]):    
            
            if starttime:
                self.data[n] = self.data[n].loc[self.data[n].index>=starttime]                
                if self.rsq:
                    self.rsq[n] = self.rsq[n].loc[self.rsq[n].index>=starttime]
                if self.NRB:
                    self.NRB[n] = self.NRB[n].loc[self.NRB[n].index>=starttime]                                       
            if endtime:
                self.data[n] = self.data[n].loc[self.data[n].index<=endtime]
                if self.rsq:
                    self.rsq[n] = self.rsq[n].loc[self.rsq[n].index<=endtime]
                if self.NRB:
                    self.NRB[n] = self.NRB[n].loc[self.NRB[n].index<=endtime]         
            if timestep:
                self.data[n] = self.data[n].resample(timestep, how = datamethod)
                if self.rsq:
                    self.rsq[n] = self.rsq[n].resample(timestep, how = datamethod)
                if self.NRB:
                    self.NRB[n] = self.NRB[n].resample(timestep, how = datamethod)
    
        if self.depolrat:
            if starttime:
                self.depolrat[0] = self.depolrat[0].loc[self.depolrat[0].index>=starttime]                                                    
            if endtime:
                self.depolrat[0] = self.depolrat[0].loc[self.depolrat[0].index<=endtime]      
            if timestep:
                self.depolrat[0] = self.depolrat[0].resample(timestep, how = datamethod)
            
        if verbose:
            print '... Done!'
        
        
        return self
        
    def range_cor(self):
        import numpy as np
        from copy import deepcopy
        
        dataout = deepcopy(self.data)        
        bg = [self.header['bg_avg2'],self.header['bg_avg1']]
      
        for n in range(self.header['numchans'][0]): 
            rsq = (np.array(dataout[n].columns, dtype=float)/1000)**2
            for i in self.data[n].index:
                dataout[n].ix[i] = (self.data[n].ix[i] - bg[n].ix[i])*rsq
        
        self.rsq = dataout        
        return self
        
    def calculate_NRB(self, showplots = False):
        
        """
        Extracts data from MPL calibration files and applies it
        to a range-corrected set of mini-MPL data to convert from counts to attenuated backscatter
        """
        import numpy as np
        import array,struct
        import os,sys
        from copy import deepcopy
        import matplotlib.pyplot as plt
        import datetime
        
        olddir = os.getcwd()
        if sys.platform == 'win32':
            newdir = 'C:\Users\dashamstyr\Dropbox\Lidar Files\MPL Data\Calibration File Archive'
        else:
            newdir = '/data/lv1/pcottle/MPLCalibration'
        os.chdir(newdir)
        
        #if data were collected before June 2013, they were collected with MPL5008 and require
        #the associated calibration files
        if self.header.index[0] < datetime.datetime(2013,6,1):
    
            deadtimefile = 'MMPL5008_deadtime.bin'
            overlapfile = 'MMPL5008_overlap.bin'
            afterpulsefile = 'MMPL5008_afterpulse.bin'
                
        
            with open(deadtimefile,'rb') as binfile:
                deadtimedat = array.array('d')
                while True:        
                    try:
                        deadtimedat.fromfile(binfile, 1)
                    except EOFError:
                        break
                
            coeffs = np.array(deadtimedat[::-1])
            
            MPLout = deepcopy(self.data)
            
            if showplots:
                temp = self.data[0].values
                maxval = temp.max()
                numsteps = maxval/100.
                x = np.arange(0,maxval,numsteps)
                y = np.polynomial.polynomial.polyval(x, coeffs)
                fig = plt.figure()
                fig.clf()
                ax = fig.add_subplot(111)
                ax.plot(x,y)
                ax.set_title('Deadtime Correction Curve')
                plt.show()
            
            
            deadtimecor = np.empty([self.header['numchans'][0],len(MPLout[0].index),len(MPLout[0].columns)])
            for n in range(self.header['numchans'][0]):
                for i in range(len(MPLout[n].index)):
                    deadtimecor[n,i,:] = np.polynomial.polynomial.polyval(MPLout[n].iloc[i],coeffs)
        
            with open(afterpulsefile, 'rb') as binfile:
                afterpulsedat = array.array('d')
                
                filedat = os.stat(afterpulsefile)
                numvals = filedat.st_size/8
                afterpulsedat.fromfile(binfile,numvals)
            
            numpairs = (numvals-1)/2
            mean_energy = np.array(afterpulsedat[0])
            aprange = np.array(afterpulsedat[1:numpairs+1])*1000.0
            apvals = np.array(afterpulsedat[numpairs+1:])
            
            if showplots:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(aprange[:100],apvals[:100])
                ax.set_title('Afterpulse Correction')
                plt.show()
            
            bg = [self.header['bg_avg2'],self.header['bg_avg1']]
              
            for n in range(self.header['numchans'][0]):
                altvals = np.array(MPLout[n].columns, dtype='float')
                interp_afterpulse = np.interp(altvals,aprange,apvals)
                
                for i in range(len(MPLout[n].index)):
                    MPLout[n].iloc[i] = MPLout[n].iloc[i]*deadtimecor[n,i] - interp_afterpulse - bg[n][i]
                    
            with open(overlapfile, 'rb') as binfile:
                overlapdat = array.array('d')
                 
                filedat = os.stat(overlapfile)
                numvals = filedat.st_size/8
                overlapdat.fromfile(binfile,numvals)
            
            numpairs = numvals/2
            overrange = np.array(overlapdat[:numpairs])*1000
            overvals = np.array(overlapdat[numpairs:])    
        
            if showplots:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(overrange,overvals)
                ax.set_title('Overlap Correction')
                plt.show()
                
            energy = self.header['energy']
            
            altvals = np.array(MPLout[0].columns, dtype='float')
            interp_overlap = np.interp(altvals,overrange,overvals)
                
            for v in range(len(altvals)):
                if altvals[v] > max(overrange):
                    interp_overlap[v] = 1.0
            
            for n in range(self.header['numchans'][0]):     
                rsq = (np.array(MPLout[n].columns, dtype=float)/1000)**2
                for i in range(len(MPLout[n].index)):
                    MPLout[n].iloc[i] = (MPLout[n].iloc[i]/(interp_overlap*energy[i]))*rsq
            
            os.chdir(olddir)
            self.NRB = MPLout
            return self
        else:        
            deadtimefile = 'MMPL5012_SPCM22625_deadtime7.bin'
            overlapfile = 'MMPL5012_Overlap_201307310000.bin'
            afterpulsefile = 'MMPL5012_Afterpulse_201308051500.mpl.bin'
    
            with open(deadtimefile,'rb') as binfile:
                deadtimedat = array.array('f')
                while True:        
                    try:
                        deadtimedat.fromfile(binfile, 1)
                    except EOFError:
                        break
                
            coeffs = np.array(deadtimedat[::-1])
            
            MPLout = deepcopy(self.data)
            
            if showplots:
                temp = self.data[0].values
                maxval = temp.max()
                numsteps = maxval/100.
                x = np.arange(0,maxval,numsteps)
                y = np.polynomial.polynomial.polyval(x, coeffs)
                fig = plt.figure()
                fig.clf()
                ax = fig.add_subplot(111)
                ax.plot(x,y)
                ax.set_title('Deadtime Correction Curve')
                plt.show()
            
            
            deadtimecor = np.empty([self.header['numchans'][0],len(MPLout[0].index),len(MPLout[0].columns)])
            for n in range(self.header['numchans'][0]):
                for i in range(len(MPLout[n].index)):
                    deadtimecor[n,i,:] = np.polynomial.polynomial.polyval(MPLout[n].iloc[i],coeffs)
        
            with open(afterpulsefile, 'rb') as binfile:
                aprange = array.array('d')  #array of range values
                apvals_copol = array.array('d')  #array of correction values for copol
                apvals_crosspol = array.array('d')  #array of correction values for crosspol
                
                
                
                temp = binfile.read(4)
                apheader = temp
                temp = binfile.read(2)
                apversion = struct.unpack('H',temp)[0] #version number for correction file
                temp = binfile.read(1)
                apnumchan = struct.unpack('B',temp)[0] #number of channels
                temp = binfile.read(4)
                apnumbins = struct.unpack('I',temp)[0] #number of correction bins
                temp = binfile.read(8)
                apenergy = struct.unpack('d',temp)[0] #laser energy during calibration in uJ
                
                if apversion != 3:
                    print 'WARNING!  This version of calibration file might be incompatible!'
                
                   
                if apnumchan == 2:
                    temp = binfile.read(8)
                    apbg_copol = struct.unpack('d',temp)[0]
                    temp = binfile.read(8)
                    apbg_crosspol = struct.unpack('d',temp)[0]
                    
                    aprange.fromfile(binfile,apnumbins)
                    aprange = [n*1000 for n in aprange] # convert range values to meters
                    apvals_copol.fromfile(binfile,apnumbins)
                    copol_norm = [n/apenergy for n in apvals_copol] #convert to counts/us
                    apvals_crosspol.fromfile(binfile,apnumbins)
                    crosspol_norm = [n/apenergy for n in apvals_crosspol] #convert to counts/us
                    apbg = [apbg_copol,apbg_crosspol]
                    apvals = [copol_norm,crosspol_norm]
                else:
                    print "Warning - wrong number of channels detected!"
            
            if showplots:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(aprange[:100],apvals[0][:100],aprange[:100],apvals[1][:100])
                ax.set_title('Afterpulse Correction')
                plt.show()
            
            bg = [self.header['bg_avg2'],self.header['bg_avg1']]
            energy= self.header['energy']
            
            for n in range(apnumchan):
                altvals = np.array(MPLout[n].columns, dtype='float') #data altitudes in m
                interp_afterpulse = np.interp(altvals,aprange,apvals[n])
                
                for i in range(len(MPLout[n].index)):
                    MPLout[n].iloc[i] = (MPLout[n].iloc[i]*deadtimecor[n,i] - interp_afterpulse - bg[n][i])/energy[i]
                    
            with open(overlapfile, 'rb') as binfile:
                overlapdat = array.array('d')
                 
                filedat = os.stat(overlapfile)
                numvals = filedat.st_size/8
                overlapdat.fromfile(binfile,numvals)
            
            numpairs = numvals/2
            overrange = np.array(overlapdat[:numpairs])*1000
            overvals = np.array(overlapdat[numpairs:])    
        
            if showplots:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(overrange,overvals)
                ax.set_title('Overlap Correction')
                plt.show()
                
            
            altvals = np.array(MPLout[0].columns, dtype='float')
            interp_overlap = np.interp(altvals,overrange,overvals)
                
            for v in range(len(altvals)):
                if altvals[v] > max(overrange):
                    interp_overlap[v] = 1.0
            
            for n in range(self.header['numchans'][0]):     
                rsq = (np.array(MPLout[n].columns, dtype=float)/1000)**2  #range in km for rsquared correction
                for i in range(len(MPLout[n].index)):
                    MPLout[n].iloc[i] = MPLout[n].iloc[i]*rsq/interp_overlap            
            
            os.chdir(olddir)
            self.NRB = MPLout
            return self
    
    def calculate_depolrat(self):
        import pandas as pan        
        copol = self.NRB[0]
        crosspol = self.NRB[1]
        
        depolMPL = crosspol.values/copol.values
        
        depolvals = depolMPL/(depolMPL+1)
        self.depolrat = []
        self.depolrat.append(pan.DataFrame(depolvals,index = copol.index, columns = copol.columns))
        
        return self
    
    def calculate_SNR(self,bg_alt=[],numprofs=1,verbose=False, datatypes=['all']):
        import pandas as pan
        import numpy as np
        """
        Calculates signal to noise ratios for mpl data
        
        inputs:
        dfin = a pandas dataframe
        bg_alt = altitude above which signal is assumed to be purely background
                 if empty, topmost 100 data points are used
        num profs = number of vertical profiles to average together, defaults to 1
        datatypes = list of data types to callculate SNR for.  Could be 
                    'raw','rsq','NRB','depolrat', or 'all'
        
        output:
        self.SNR = a dict of pandas dataframes with datatype keys containing 
                    SNR values
        
        """
        bg = [self.header['bg_avg2'],self.header['bg_avg1']]
        
        if verbose:
            print "Calculating SNR"
        
        datasets=[]
        SNRdict={}
        for d in datatypes:
            if d=='data' or d=='all':
                datasets.append(('data',self.data))
            if d=='rsq' or d=='all':
                if self.rsq:
                    datasets.append(('rsq',self.rsq))
                elif verbose:
                    print "No RSQ data available for SNR calc"
            if d=='NRB' or d=='all':
                if self.NRB:
                    datasets.append(('NRB',self.NRB))
                elif verbose:
                    print "No NRB data available for SNR calc"
            if d=='depolrat' or d=='all':
                if self.depolrat:
                    datasets.append(('depolrat',self.depolrat))
                elif verbose:
                    print "No depolrat available for SNR calc"
        
        for dset_name,dset in datasets:            
            for n in range(len(dset)): 
                tempdat=dset[n]
                if not bg_alt:
                    bg_alt=tempdat.columns[-100]
                SNRtemp=pan.DataFrame(np.empty_like(tempdat.values),index=tempdat.index,columns=tempdat.columns)
                for i in tempdat.index:
                    tempprof=tempdat.ix[i]
                    if dset_name=='data':
                        tempback=bg[n].ix[i]
                    else:
                        tempback=0
                    tempvals=[v for v,r in zip(tempprof.values,tempprof.index) if r>=bg_alt]
                    tempfilt=[x for x in tempvals if not np.isnan(x)]
                    sigmatemp=np.std(tempfilt)
                    Ctemp=sigmatemp/np.mean(np.sqrt(np.abs(tempfilt)))
                    SNR = lambda x: (x-tempback)/(Ctemp*np.sqrt(np.abs(x)))
                    
                    SNRprof=np.array([SNR(v) for v in tempprof.values]).clip(0)
                    SNRtemp.ix[i]=SNRprof
                SNRdict[dset_name][n]=SNRtemp
                
        self.SNR=SNRdict

        if verbose:
            print "SNR calculation done!"
            
        return self
    
    def calc_all(self):
        """
        calculates all uncalculated fields for an mpl object
        """
        if not self.rsq:
            self.range_cor()
        
        if not self.NRB:
            self.calculate_NRB()
        
        if not self.depolrat:
            self.calculate_depolrat()
        
        if not self.SNR:
            self.calculate_SNR()
        
        return self
    
           
      

if __name__ == '__main__':
    
    import os
    
    olddir = os.getcwd()
    
    os.chdir('C:\\SigmaMPL\DATA')
    
    filename = '201312010000.mpl'
    
    print 'Testing MPLtoHDF'
    MPLtoHDF(filename)    
    print 'done'
    
    print 'Testing MPL class functions'
    
    print 'Import MPL data from .mpl file'
    
    MPLtest = mtools.MPL()
    MPLtest.fromMPL(filename)
    
    print 'Done'
    
    print 'Resample in altitude and time'
#    
    
    
    os.chdir(olddir)
