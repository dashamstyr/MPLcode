import pandas as pan
import numpy as np
import os
import MPLtools as mtools
import image_process_tools as imtools
import matplotlib.pyplot as plt
import MPL_plot as mplot

#----------------------------------------------------------------------------
#Uses tools created in MPL_tools to open all files in a folder and resample
#them to a regular spacing in altitude/date the concatenates them into one
#pandas dataframe and plots it using LNC_plot
#July 05, 2012
#----------------------------------------------------------------------------

os.chdir('C:\SigmaMPL\DATA')

#newdir = mtools.set_dir('Select Event Folder')
#
#os.chdir(newdir)
#
#files = os.listdir(newdir)

files = mtools.get_files('Select MPL files', filetype = ('.mpl', '*.mpl'))
rawfiles = []

#set altitude range and date step sizes

altrange = np.arange(150,20000,10)#meters
timestep = '60S' #seconds

#check to see if each file has been processed before and separate processed
#files into a new list

for f in files:
    if '.mpl' in f:
        rawfiles.append(f)

#open, altitude resample, and concatenate data and mask files

print rawfiles
for r in rawfiles:
    MPLdat_temp = mtools.MPL()
    MPLdat_temp.fromMPL(r)
    MPLdat_temp.alt_resample(altrange)

    try:
        MPLdat_event.append(MPLdat_temp)
    except NameError:
        MPLdat_event = MPLdat_temp
   
#sort by index to make certain data is in order then set date ranges to match
for n in range(MPLdat_event.header.numchans):
    data = MPLdat_event.data[n]
    data = data.sort_index()

start = data.index[0]
end = data.index[-1]

MPLdat_event.time_resample(timestep, timerange = [start,end])

MPLdat_event.background_subtract()

MPLdat_event.range_cor()

MPLdat_event.power_cor()


#for n in range(MPLdat_event.header.numchans):
#    MPLdat_event.data[n] = imtools.blur_image(MPLdat_event.data[n],(5,5), kernel ='Flat')

olddir = os.getcwd()

try:
    os.chdir('Processed')
except WindowsError:
    os.makedirs('Processed')
    os.chdir('Processed')

if len(rawfiles) == 1:   
    [path,startfile] = os.path.split(rawfiles[0])
    d_filename = startfile.split('.')[0]+'_proc'
else:
    [path,startfile] = os.path.split(rawfiles[0])
    [path,endfile] = os.path.split(rawfiles[-1])
    d_filename = startfile.split('.')[0]+'-'+endfile.split('.')[0]+'_proc'

print 'Saving '+d_filename
MPLdat_event.save_to_HDF(d_filename+'.h5')
#MPLdat_event.save_to_HDF(d_filename+'.h5')
print 'Done'
os.chdir(olddir)

try:
    os.chdir('Figures')
except WindowsError:
    os.makedirs('Figures')
    os.chdir('Figures')

#create figure and plot image of depolarization ratios
fsize = 18 #baseline font size
ar = 2.0  #aspect ratio
figheight = 12 #inches

plt.rc('font', family='serif', size=fsize)


fig = plt.figure()
fig.clf()
h_set = range(1,25)
h_set = map(str,h_set)

copol = MPLdat_event.data[0]
crosspol = MPLdat_event.data[1]

depolMPL = crosspol.values/copol.values

depolvals = depolMPL/(depolMPL+1)

depol = pan.DataFrame(depolvals,index = copol.index, columns = copol.columns)
#depol = imtools.blur_image(depol,(7,7), kernel ='Flat')

datetime = copol.index
alt = copol.columns

print 'Generating Figure'

ax1 = fig.add_subplot(2,1,1)
im1 = mplot.backscatter_plot(fig, ax1, ar, datetime,alt[::-1],copol.T[::-1], (0,0.003), fsize = fsize)
cbar1 = fig.colorbar(im1, orientation = 'vertical', aspect = 6)
cbar1.ax.tick_params(labelsize = fsize)
mplot.dateticks(ax1, datetime, hours = h_set, fsize = fsize, tcolor = 'w')
ax1.set_xticklabels([])
t1 = ax1.set_title('Attenuated Backscatter', fontsize = fsize+10)
t1.set_y(1.03)
        
ax2 = fig.add_subplot(2,1,2)
im2 = mplot.depol_plot(fig, ax2, ar, datetime,alt[::-1],depol.T[::-1], (0,0.5), fsize = fsize)
cbar2 = fig.colorbar(im2, orientation = 'vertical', aspect = 6)
cbar2.set_ticks(np.arange(0,0.6,0.1))
cbar2.set_ticklabels(np.arange(0,0.6,0.1))
cbar2.ax.tick_params(labelsize = fsize)
#set axis ranges and tickmarks based on data ranges
mplot.dateticks(ax2, datetime, hours = h_set, fsize = fsize)
ax2.set_xlabel('Time [Local]',fontsize = fsize+4)
fig.autofmt_xdate()
t2 = ax2.set_title('Linear Depolarization Ratio',fontsize = fsize+10)
t2.set_y(1.03)    

##plt.savefig(savetitle,dpi = 100, edgecolor = 'b', bbox_inches = 'tight')
fig.set_size_inches(figheight*ar,figheight) 
plt.savefig(d_filename+'.png')
print 'Done'


os.chdir(olddir)