'''
Created on Apr 20, 2012

@author: jeven
'''
import TableIO
import numpy as np
import os
import errno
import time

from UtilityFunctions import *
from CAMB_Config import *

def ImportTransferFunction(transfer_file):
    """
    Imports the Transfer Function file to be analysed, and returns the pair ln(k), ln(T)
    
    Input: "transfer_file": full path to the file containing the transfer function (from camb).
    
    Output: ln(k), ln(T)
    """
        
    print "  Importing Transfer Function"
    k,T = TableIO.readColumns(transfer_file,"!#",columns=[0,1])
        
    k = np.log(k)
    T = np.log(T)
                
    if config.loud:
        print "    Minimum logk "+str( min(k))
        print "    Maximum logk "+str( max(k))

    return k,T
        
def MassSetup(min_M,max_M,M_step):
    """
    Simply initializes a vector of radius values from input data
    
    Input: None (all necessary inputs defined in config)
    Output: R: an array of values corresponding to different values of the radius.
    """
        
    M = np.arange(min_M,max_M,M_step)
    return M


def InputDirectorySetup():
    """
    Checks the directory structure for input.
    
    Input: None (all relevant parameters defined in config)
    Output: transfer_file: full path to the transfer function file.
    """

    i=0
    while not os.path.isdir(config.camb_folder) and i<4:
        print "Your cam folder is currently: "+config.camb_folder
        print "Error: the CAMB folder specified in config.py is incorrect. "
        i = i+1
        config.camb_folder = input("Please specify it now:")
            
    if i==4:
        sys.exit("You failed to enter the correct directory for CAMB, please try in the config file.")
            
    SetParameter(os.path.abspath(os.path.dirname(sys.argv[0]))+"/Config.py",{"camb_folder":config.camb_folder})
        
    if not os.path.exists(config.camb_folder+"/camb"):
        sys.exit("Your CAMB is not installed in the usual way, please reinstall it")
            
    if not os.path.exists(config.camb_folder+'/'+config.transfer_file):
        print "No power spectrum file exists under your CAMB folder with the given filename:"
        print "Using CAMB to produce one (CURRENTLY ONLY WITH DEFAULT/CURRENT PARAMS)"
        CAMB()
        config.use_camb = False
            
                    
    transfer_file = config.camb_folder+'/'+config.transfer_file              
    if not os.path.exists(transfer_file):
        sys.exit("  Could not find transfer function file. Please enter it correctly in config.")   
         
    return transfer_file
  
def CAMB():
    """
    Uses CAMB in its current setup to produce a transfer function
    
    The function needs to be imported by calling ImportTransferFunction.
    """

    os.chdir(config.camb_folder)

    SetParameter(config.camb_folder+'/params.ini', camb_config.camb_dict)
    os.system('./camb params.ini')
        
    config.transfer_file = ReadParameter( config.camb_folder+'/params.ini', "output_root")+"_transfer_out.dat"    
    SetParameter( os.path.abspath(os.path.dirname(sys.argv[0]))+"/Config.py", {"transfer_file":config.transfer_file})
        
        
def OutputDirectorySetup():
    """
    Sets up the output directories.
    
    Input: none (all defined in config)
    Output: output_prefix: a directory labelled by the current time inside the project directory (config)
    """

    i=0
    while not os.path.isdir(config.project_dir) and i<4:
        print "No directory "+config.project_dir+" found in this location"
        location = input("Would you like to create it (0) or specify an existing/new project location (1)?")
        if location == 0:
            try:
                os.makedirs(config.project_dir)
            except OSError, e:
                if e.errno != errno.EEXIST:
                    raise
        elif location == 1:
            config.project_dir = input("Please enter the location for your project: ")
            try:
                os.makedirs(config.project_dir)
            except OSError, e:
                if e.errno != errno.EEXIST:
                    raise
        else:
            print "Please enter either 0 or 1:"
        
        i = i+1
        
    if config.create_time_subdirs:
        output_prefix = config.project_dir+'/'+time.strftime("%Y-%m-%d--%H%M")
        try:
            os.makedirs(output_prefix)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
    else:
        output_prefix = config.project_dir
    
    if config.output_table:
        try:
            os.makedirs(output_prefix+'/OutputTables')
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
            
    config.out_pref = output_prefix        
    SetParameter(os.path.abspath(os.path.dirname(sys.argv[0]))+"/Config.py", {"project_dir":config.project_dir})   
       

def Setup():
    """
    A convenience function used to fully setup the workspace in the 'usual' way
    """
    
    transfer_file = InputDirectorySetup()
    if config.use_camb:
        CAMB()
        
    OutputDirectorySetup()
    M = MassSetup(config.min_M,config.max_M,config.M_step)
    k,T = ImportTransferFunction(transfer_file)
    
    return k,T,M
    