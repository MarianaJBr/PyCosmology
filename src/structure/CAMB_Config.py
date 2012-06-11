'''
Created on Apr 13, 2012

@author: smurray
'''

class camb_config(object):
    """
    Defines some paramaters to be used in the CAMB run (if used). 
    
    The reason they are available here, rather than just setting them in the params.ini file of the camb folder is twofold.
    Firstly, only 'important' parameters for the power spectrum modelling are included here, and also the labelling and information
    given is cleaner and easier to see. Together this makes this file easier to change for power spectrum/mass variance purposes.
    Secondly, in the future if parameter space is to be searched, it will be easier modifying the parameters through this program.
    """
    #=====================================================================================
    # LINEAR OR NON-LINEAR POWER SPECTRUM?
    #=====================================================================================
    # 0: linear, 1: non-linear matter power (HALOFIT), 2: non-linear CMB lensing (HALOFIT)
    do_nonlinear = 0
    
    #=====================================================================================
    # MAIN COSMOLOGICAL PARAMETERS
    #=====================================================================================
    use_physical   = False
    
    if use_physical:
        ombh2          = 0.03
        omch2          = 0.112
        omnuh2         = 0
        omk            = 0
    else:
        omega_baryon   = 0.0456
        omega_cdm      = 0.2274
        omega_lambda   = 0.727
        omega_neutrino = 0
    hubble         = 70.5
    # Effective equation of state parameter for dark energy, assumed constant
    w              = -1
    
    # Constant comoving sound speed of the dark energy (1=quintessence)
    cs2_lam        = 0
     
    temp_cmb           = 2.725
    helium_fraction    = 0.24

    #=====================================================================================
    # INITIAL CONDITIONS
    #=====================================================================================
    #Initial scalar perturbation mode (adiabatic=1, CDM iso=2, Baryon iso=3, 
    # neutrino density iso =4, neutrino velocity iso = 5) 
    initial_condition   = 1


    #=====================================================================================
    # TRANSFER/POWER SPECTRUM DETAILS
    #=====================================================================================    
    #Transfer function settings, transfer_kmax=0.5 is enough for sigma_8
    #transfer_k_per_logint=0 sets sensible non-even sampling; 
    #transfer_k_per_logint=5 samples fixed spacing in log-k
    #transfer_interp_matterpower =T produces matter power in regular interpolated grid in log k; 
    # use transfer_interp_matterpower =F to output calculated values (e.g. for later interpolation)
    transfer_high_precision = True
    transfer_kmax           = 1000
    transfer_k_per_logint   = 50
    transfer_num_redshifts  = 1
    transfer_interp_matterpower = True

    #=====================================================================================
    # ACCURACY/PERFORMANCE OPTIONS
    #=====================================================================================
    #Default scalar accuracy is about 0.3% (except lensed BB) if high_accuracy_default=F
    #If high_accuracy_default=T the default taget accuracy is 0.1% at L>600 (with boost parameter=1 below)
    #Try accuracy_boost=2, l_accuracy_boost=2 if you want to check stability/even higher accuracy
    #Note increasing accuracy_boost parameters is very inefficient if you want higher accuracy,
    #but high_accuracy_default is efficient 
    
    high_accuracy_default=True
    
    #Increase accuracy_boost to decrease time steps, use more k values,  etc.
    #Decrease to speed up at cost of worse accuracy. Suggest 0.8 to 3.
    accuracy_boost          = 1
    
    #Larger to keep more terms in the hierarchy evolution. 
    l_accuracy_boost        = 1
    
    #Increase to use more C_l values for interpolation.
    #Increasing a bit will improve the polarization accuracy at l up to 200 -
    #interpolation errors may be up to 3%
    #Decrease to speed up non-flat models a bit
    l_sample_boost          = 1
    
#*****************************************************************************************
#*****************************************************************************************
# DON'T CHANGE THIS BELOW    
#*****************************************************************************************
#*****************************************************************************************
   
    if use_physical:
            
        camb_dict = {"do_nonlinear" : do_nonlinear,
                    "use_physical"   : use_physical,
                    "w"              : w,
                    "cs2_lam"        : cs2_lam,
                    "temp_cmb"       : temp_cmb,
                    "helium_fraction": helium_fraction,
                    "initial_condition": initial_condition,
                    "transfer_high_precision":transfer_high_precision,
                    "transfer_kmax"  :transfer_kmax,
                    "transfer_k_per_logint":transfer_k_per_logint,
                    "transfer_num_redshifts"  : transfer_num_redshifts,
                    "transfer_interp_matterpower" : transfer_interp_matterpower,
                    "high_accuracy_default":high_accuracy_default,
                    "accuracy_boost"  : accuracy_boost,
                    "l_accuracy_boost": l_accuracy_boost,
                    "l_sample_boost"  : l_sample_boost,
                    "ombh2"          : ombh2,
                    "omch2"          : omch2,
                    "omnuh2"         : omnuh2,
                    "omk"            : omk,
                    "hubble"         : hubble}
            
    else:
            
        camb_dict = {"do_nonlinear" : do_nonlinear,
                    "use_physical"   : use_physical,
                    "w"              : w,
                    "cs2_lam"        : cs2_lam,
                    "temp_cmb"       : temp_cmb,
                    "helium_fraction": helium_fraction,
                    "initial_condition": initial_condition,
                    "transfer_high_precision":transfer_high_precision,
                    "transfer_kmax"  :transfer_kmax,
                    "transfer_k_per_logint":transfer_k_per_logint,
                    "transfer_num_redshifts"  : transfer_num_redshifts,
                    "transfer_interp_matterpower" : transfer_interp_matterpower,
                    "high_accuracy_default":high_accuracy_default,
                    "accuracy_boost"  : accuracy_boost,
                    "l_accuracy_boost": l_accuracy_boost,
                    "l_sample_boost"  : l_sample_boost,
                    "omega_baryon"   : omega_baryon,
                    "omega_cdm"      : omega_cdm,
                    "omega_lambda"   : omega_lambda,
                    "omega_neutrino" : omega_neutrino}