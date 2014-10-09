# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 14:01:06 2014

@author: User
"""

import pandas as pan
from Tkinter import Tk
import tkFileDialog
import re
import pandas as pan

def get_files(titlestring,filetype = ('.h5','*.h5')):
     
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
   
    root.destroy()
    return result

def h5extractor(filename=[]):  
    if not filename:
        filename=get_files('Select one HDF5 file to be parsed')
        copoldat_NRB = pan.read_hdf(filename[0],'copol_NRB')
    else:
        copoldat_NRB = pan.read_hdf(filename,'copol_NRB')
    
    return copoldat_NRB
    
if __name__=='__main__':
    x=h5extractor()