'''
Created on Apr 20, 2012

@author: jeven
'''

if __name__ == '__main__':


#===================================================================
#    IMPORT LOCAL MODULES AND CONFIG FILES
#===================================================================
    import MassFunctionClasses as mfc
    import UtilityFunctions as ute
    import SetupFunctions as setup
    import OutputFunctions as output
    from Observation import *
    from Config import *

#===================================================================
#    IMPORT EXTERNAL MODULES
#===================================================================    
    import time
    import numpy as np
    
    total_time_0 =time.time()
#===================================================================
#    SETUP THE CALCULATIONS
#===================================================================     
    # Write important parameter values to stdout used in calculation.
    ute.WriteHead()
    
    print "SETTING UP"
    sec_time_0 = time.time()
    # Use the convenience function Setup() to import the transfer
    # function, set up a vector of radii, establish a step-size for
    # Romberg integration, and create output directories.
    # If use_camb is true, then CAMB operations will also be performed
    k_vector, transfer, masses = setup.Setup()

    masses = 10**masses
    sec_time_1 = time.time()
    print "  Section Time: ", sec_time_1-sec_time_0
#===================================================================
#    FIND RELEVANT QUANTITIES
#=================================================================== 
    # Initialize the perturbation calculations with the transfer
    # function, normalizing and extrapolating (or truncating).
    print "PERFORMING CALCULATIONS"
    sec_time_0 = time.time()
    mass_var = np.zeros((len(config.z),len(config.WDM)+1,masses.shape[0]))
    mass_func = np.zeros((len(config.approach),len(config.z),len(config.WDM)+1,masses.shape[0]))
    radii = np.zeros_like(masses)
    if config.output_dlnsdlnm:
        dlnsdlnm = np.zeros((len(config.z),len(config.WDM)+1,masses.shape[0]))
    if config.output_dndm:
        dndm = np.zeros((len(config.approach),len(config.z),len(config.WDM)+1,masses.shape[0]))
    if config.output_nufnu:
        nufnu = np.zeros((len(config.z),len(config.WDM)+1,masses.shape[0]))
    
    
    for i,wdm in enumerate([None]+config.WDM):
        for j,z in enumerate(config.z):
            print "  Using Cosmology with WDM Mass", wdm, "keV, z = ", z
                
            pert = mfc.Perturbations(k_vector,transfer,z,config.step_size,config.cosmology,wdm,[config.k_begins_at,config.k_ends_at])
            if j==0:
                if i == 0:
                    power_spec = np.atleast_2d(pert.power_spectrum)
                else:
                    power_spec = np.append(power_spec,[pert.power_spectrum],axis=0)
         
            for k,approach in enumerate(config.approach): 
                print "  And the", approach, "Approach"
                for ii, M in enumerate(masses):     
                    
     
                    # Find various quantities via perturbation class for vector of R.  
                    if config.output_radius:
                        radii[ii] = pert.Radius(M)
                    mass_var[j,i,ii] = pert.MassVariance(M)
                    
                    if config.output_dlnsdlnm:
                        dlnsdlnm[j,i,ii] = pert.dlnsdlnM(M)
                    if config.output_dndm:
                        dndm[k,j,i,ii] = pert.dndM(M,approach)
                    if config.output_nufnu:
                        nufnu[j,i,ii] = pert.nufnu(M,approach)                    
                        

                    mass_func[k,j,i,ii] = pert.MassFunction(M,approach)  
                    
                    
                      
                if config.output_table:
                    output_data = [] 
                    output_data = output_data + [np.log10(masses)]   
                    if config.output_radius:
                        output_data= output_data+[radii]
                
                    if config.output_massvar:
                        output_data= output_data+[np.sqrt(mass_var)[j,i,:]]
                
                    if config.output_dlnsdlnm:
                        output_data= output_data+[dlnsdlnm[j,i,:]]
                        
                    if config.output_nufnu:
                        output_data= output_data+[nufnu[j,i,:]]
                        
                    if config.output_dndm:
                        output_data= output_data+[dndm[k,j,i,:]]
                        
                    if config.output_ps:
                        output_data= output_data+[mass_func[k,j,i,:]]
                     
                    output.EverythingOutput(output_data,z,approach,wdm = wdm)
         
    sec_time_1 = time.time()
    print "  Section Time: ", sec_time_1-sec_time_0   
    
#============================================================
#    SIMULATION COMPARISON
#============================================================
    hist = []
    mass_bins = []
    if config.compare_sims:
        for i,filename in enumerate(config.sim_files):
            obs = HMFObservations(config.sim_dir+filename,20)
            (hist_temp,mass_bins_temp) = obs.MakeHist(config.hist_steps)
            mass_bins_temp_2 = np.zeros((mass_bins_temp.shape[0]-1))
            for j in range(mass_bins_temp.shape[0]-1):
                mass_bins_temp_2[j] = mass_bins_temp[j]+(mass_bins_temp[j+1] - mass_bins_temp[j])/2.0
            hist = hist + [hist_temp]
            mass_bins = mass_bins+[mass_bins_temp_2]

    
#============================================================
#    PLOT RESULTS
#============================================================
    print "PLOTTING RESULTS"
    sec_time_0 = time.time()
    if config.plot_power:
        output.PowerPlot(pert.k_extended,power_spec)
        
    if config.plot_massvar:
        output.SigmaPlot(np.log10(masses), mass_var)
        output.SigmaLogPlot(np.log10(masses), mass_var)
        
    if config.plot_massfunc:
        if not config.compare_sims:
            output.MassFunctionPlot(np.log10(masses), mass_func,sims=None,mass_bins = None)
        else:
            output.MassFunctionPlot(np.log10(masses), mass_func,sims=hist,mass_bins=mass_bins)

    
    sec_time_1 = time.time()
    print "  Section Time: ", sec_time_1-sec_time_0  
        
    total_time_1 = time.time()
    print "Total Time: ", total_time_1 - total_time_0    

            
        
    
    
    
    