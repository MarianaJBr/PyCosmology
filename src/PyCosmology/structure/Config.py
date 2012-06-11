class config(object):
    """
    Defines parameters to be used in main.py
    """
    
#===================================================================
#    DIRECTORY SETUP
#===================================================================     
    # Define the name of the input matter power spectrum file (from CAMB)
    transfer_file = 'test_transfer_out.dat'
    
    # Setup Directories for output/input
    project_dir = '/Users/smurray/Documents/PhD/OverallProject/MassVarianceOutput'
    
    #If true, all file output will be under a folder labeled by the time of the run.
    # I hope to later add an option to label the folder by parameters.
    create_time_subdirs = True 
                               
    camb_folder = '/Users/smurray/Documents/PhD/OverallProject/ExternalCode/camb'

 
#===================================================================
#    ACCURACY PARAMETERS
#===================================================================   
    # Define the step size for the integrated function in log space.
    step_size = 0.015625
 
    # Extrapolation parameters.
    extrapolate = True
    k_ends_at = 2000.0 #Lower cut-off corresponds to small-scale damping.
    k_begins_at = 0.00000001  # Higher cut-off corresponds to finite volume damping in simulations.   
    
#===================================================================
#    RUN PARAMETERS
#===================================================================     
    # Which values of the radius to use?
    min_M = 8.0
    max_M = 16.0
    M_step = 0.15
 

    # Whether to run CAMB and use its transfer function on-the-fly
    # (See CAMB_Config)   
    use_camb = False
    
    # Mass Function fit
    approach = ["Reed07"]
#===================================================================
#   TEST PARAMETERS
#===================================================================     
    # Whether to perform a test of step size effects on integration.
    step_size_analysis = False
    step_size_number = 10
    
    # Define whether to write out extra information
    loud = False
    
#===================================================================
#   COSMOLOGICAL PARAMETERS
#=================================================================== 
    # Mean Density of the Universe (by Omega*h^2*M_sun/Mpc^3)
    mean_dens = 0.273*2.7755*10**11 
    
    # Critical Overdensity corresponding to spherical collapse
    crit_dens = 1.686
    
    # Power spectral index
    n = 1.0
    
    # Mass variance on a scale of 8Mpc
    sigma_8 = 0.812
    
    # Redshift at which to calculate the mass variance.
    z = [0.0]
    
    # Hubble parameter (make sure camb config is modified also)
    hubble = .705
    
    # WDM particle masses (empty list if none)
    WDM = []
#===================================================================
#   SIMULATION COMPARISON
#=================================================================== 
    compare_sims = False
    sim_dir = "/users/smurray/Documents/PhD/OverallProject/Simulations/"
    sim_files = ["snap_100.cdm.gcat",
                 "snap_100.wdm0.5.gcat",
                 "snap_100.wdm1.gcat",
                 "snap_100.wdm2.gcat"]    
    
    hist_steps = 0.25
#===================================================================
#   OUTPUT OPTIONS
#===================================================================        
    plot_power = True
    plot_massvar = True
    plot_massfunc = True
    
    output_table = True
    output_mass = True
    output_nufnu = True
    output_radius = True
    output_massvar = True
    output_dlnsdlnm = True
    output_dndm = True
    output_ps = True


#===================================================================
#   GLOBAL VARIABLES - DON'T TOUCH
#=================================================================== 
    out_pref = ''
    
#===================================================================
#   PARAMETER COMPILATION - DON't TOUCH
#=================================================================== 
    cosmology = {"sigma_8":sigma_8,"n":n,"mean_dens":mean_dens,"crit_dens":crit_dens,"hubble":hubble}