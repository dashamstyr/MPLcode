# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 13:02:55 2014

@author: dashamstyr
"""

import os,sys

def locations(historyfile = "MPL_history_file.txt"):
    
    with open(historyfile, 'r') as hist:
        lines = hist.readlines()
        
        locations = []
        start_dates = []
        end_dates = []
        
        for l in lines:
            if l.endswith('\r\n'):
                l=l[:-2]
            if l.endswith('\n'):
                l=l[:-1]
            if "Location" in l:
                locations.append(l.split(":")[1])
            if "Start_date" in l:
                start_dates.append(l.split(":")[1])
            if "End_date" in l:
                if l.split(":")[1]:
                    end_dates.append(l.split(":")[1])
                else:
                    end_dates.append("na")
        
        headerdict = dict(zip(locations,zip(start_dates,end_dates)))
    
    return headerdict
    
def listdates(mplfiles):
    filenames=[]
    dates = []
    
    for f in sorted(mplfiles):
        filenames.append(f)
        tempyear = f[:4]
        tempmonth = f[4:6]
        tempday = f[6:8]
        dates.append((tempyear,tempmonth,tempday))
    
    return filenames,dates
                
def findloc(mplfile,locdict):
    filedate = int(mplfile[:8])
    
    for loc,dates in locdict.iteritems():
        startdate = int(dates[0])
        try:
            enddate = int(dates[1])        
        except ValueError:
            enddate = '99999999'
            
        if startdate <= filedate <= enddate:
            return loc   
        
    