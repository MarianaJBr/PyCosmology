'''
Created on Apr 20, 2012

@author: jeven
'''

import sys
from Config import *
from CAMB_Config import *

def ReadParameter(filename,parameter,value_type="string"):
    """
    Reads and sets a parameter from a config file (notably CAMB)
    """

    print "    Importing the parameter "+parameter+" from the file: "+filename
        
    file_object = open(filename,'r')
    for line in file_object.readlines():
        if line.strip().startswith(parameter):
            param =  line.partition('=')[2].strip()
            break
        
    if value_type is "string":
        return param
    elif value_type is "int":
        return int(param)
    elif value_type is "float":
        return float(param)
    elif value_type is "bool":
        return bool(param)
    else:
        sys.exit("Incorrect type input to method ReadParameter")
            
    file_object.close
        
def SetParameter(filename,keys):
    """
    Sets the relevant parameters, set by 'keys' into the config file 'filename'
    """
        
    print "  Setting Various Parameters in the File:"
    print "  "+filename
        
    file_object = open(filename,'r+')
        
    file_data = file_object.readlines()
    for parameter,value in keys.iteritems():
        for number,line in enumerate(file_data):
            if line.strip().startswith(parameter):
                if type(value) is type ("string"):
                    if filename.endswith("params.ini"):
                        file_data[number] = line.replace(line.partition('=')[2].strip(),value) 
                    else:
                        file_data[number] = line.replace(line.partition('=')[2].strip(),"'"+value+"'") 
                elif type(value) is type(True):
                    if filename.endswith("params.ini"):
                        if value:
                            file_data[number] = line.replace(line.partition('=')[2].strip(),'T')
                        else:
                            file_data[number] = line.replace(line.partition('=')[2].strip(),'F')
                    else:
                        file_data[number] = line.replace(line.partition('=')[2].strip(),str(value)) 
                else:
                    file_data[number] = line.replace(line.partition('=')[2].strip(),str(value))
                    
                break
                
                
    file_object.seek(0)
    file_object.writelines(file_data)
    file_object.close()
    
def WriteHead():
    """
    Writes the heading with appropriate information to screen and logfile
    """

    print ""
    print ""
    print "THIS PROGRAM CALCULATES THE MASS VARIANCE FOR A RANGE "
    print "--OF RADII GIVEN A MATTER POWER SPECTRUM FROM CAMB --"
    print "*****************************************************"
    print "LOCAL PARAMETERS USED IN THIS RUN"
    print ""
    print "  Approximate Step-size of the Integration: "+str(config.step_size)
    if config.extrapolate:
        print "  Power Spectrum extrapolated between "+str(config.k_begins_at)+" and "+str(config.k_ends_at)
    print ""
    print "  Loud: "+ str(config.loud)
    print ""
        
    print "EXTERNAL PARAMETERS USED IN POWER SPECTRUM ESTIMATION"
    print ""
        
    print "  do_nonlinear               :"+str(camb_config.do_nonlinear) 
    print "  use_physical               :" +str(camb_config.use_physical)
    print "  w:                         :"              +str( camb_config.w)
    print "  cs2_lam:                   :"        +str( camb_config.cs2_lam)
    print "  temp_cmb:                  :"       +str( camb_config.temp_cmb)
    print "  helium_fraction:           :"+str( camb_config.helium_fraction)
    print "  initial_condition          :"+str( camb_config.initial_condition)
    print "  transfer_high_precision    :"+str(camb_config.transfer_high_precision)
    print "  transfer_kmax              :"  +str(camb_config.transfer_kmax)
    print "  transfer_k_per_logint      :"+str(camb_config.transfer_k_per_logint)
    print "  transfer_num_redshifts     :"  +str( camb_config.transfer_num_redshifts)
    print "  transfer_interp_matterpower:" +str( camb_config.transfer_interp_matterpower)
    print "  high_accuracy_default      :"+str(camb_config.high_accuracy_default)
    print "  accuracy_boost             :"  +str( camb_config.accuracy_boost)
    print "  l_accuracy_boost           :"+str( camb_config.l_accuracy_boost)
    print "  l_sample_boost             :"  +str( camb_config.l_sample_boost)
    if camb_config.use_physical:
        print "  ombh2                      :" +str( camb_config.ombh2)
        print "  omch2                      :"  +str( camb_config.omch2)
        print "  omnuh2                     :"  +str( camb_config.omnuh2)
        print "  omk                        :"    +str(camb_config.omk)
        print "  hubble                     :"  +str(camb_config.hubble)
    else:
        print "  omega_baryon               :"+str(camb_config.omega_baryon)
        print "  omega_cdm                  :"  +str(camb_config.omega_cdm)
        print "  omega_lambda               :"  +str(camb_config.omega_lambda)
        print "  omega_neutrino             :"+str(camb_config.omega_neutrino)
    print "*****************************************************"
