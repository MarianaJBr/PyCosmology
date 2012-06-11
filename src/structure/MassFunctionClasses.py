'''
Created on Apr 5, 2012

@author: smurray
'''

#!/usr/bin/python

history = \
"""
Version 1.3    13/04/2012    Added: Writing out of logfile as well as stdout.
                                    Could be a lot cleaner (as could stdout)

Version 1.2    13/04/2012    Added: Ability to truncate the power spectrum at any point
                                    in k-space for use in either extended accuracy or 
                                    in testing low or high-k damping.
                                    
                             Added: Ability to set 'useful' parameters in CAMB via the
                                    CAMB_Config file and run them to use the gained
                                    power spectrum
                                    
                             Added: Timing for performance verifcation.
                             
Version 1.1    12/04/2012    Added: Ability to perform a step-size test to determine the
                                    optimal step-size for accuracy and speed. For the 
                                    power spectrum used, the optimal step-size turns out
                                    to be 0.0026139891392, which is used if a step_size
                                    of 0.001953125 is input to config.
                                    
Version 1.0:  10/04/2012    First working version. Currently able to calculate the mass
                            variance for a range of radii, saving a plot of the original 
                            power spectrum, and the mass variance function (logarithmic
                            and non-logarithmic). The power spectrum can be supplied (as
                            previously run CAMB output), but if it is not available,
                            the program will call CAMB and produce it on-the-fly.
                            
                            It (semi-) intelligently handles folders and working directories
                            so that the program can be run from anywhere and SHOULD work.
                            To ease the pain of supplying directories, they are are written
                            to config as defaults for all runs thereafter.
                            
                            What it does not do (and I would like it to do): be able to manage
                            (some of) the parameters that CAMB uses within the body of the 
                            program itself. Write to folders labelled by parameter values.
                            Produce nicer looking plots. Have a GUI. Test efficiency/performance
                            for different step sizes.
              
"""



from scipy.interpolate import interp1d
import scipy.integrate as intg
import numpy as np
import time
import sys

from Config import *
import SetupFunctions as setup
from CAMB_Config import *




class Perturbations(object):
    """
    A class which contains the functions necessary to evaluate the HMF and Mass Variance.
    
    Required Input:
        k: numpy vector of wavenumbers for the Transfer function
        Transfer: numpy vector of values for the Transfer function
        
    Optional Input:
        z: a float or vector of floats containing the redshift(s) of the analysis.
            Default z = 0.0
        M: a float or vector of floats containing the log10(Solar Masses) at which to perform analyses.
            Default M = 10^11 solar mass
        step_size: a float or vector of floats containing the step size(s) of the integration
                    in log space.
                    Default step_size = 0.001
        cosmology: a dictionary of cosmological parameters. Entries may include omega_lambda,
                   omega_m, sigma_8, n (initial power spectral index).
                   Default cosmology = {"sigma_8":0.812,"n":1}
        k_bounds: a list defining two values: the lower and upper limit of ln(k). Used to truncate/extend
                  the power spectrum. Default is empty, ie. no truncation or extension.
        WDM: a list of warm dark matter particle sizes in keV. Default is empty
                  
    The purpose and idea behind the class is to be able to calculate many quantities associated
    with perturbations in the early Universe. The class is initialized to form a cosmology and
    takes in various options as to how to calculate the various quantities. Most non-initial 
    methods of the class need no further input after the cosmology has been initialized.
    

    
    Non-Initial Methods:
    
    1. Top Hat Window Function: Produces a top hat window function in Fourier Space

    2. Mass Variance          : Produces the mass variance as a function of radius/mass.

    3. dlnsdlnM               : Finds the derivative of the log of sigma with respect to log mass.

    4. dndM                   : the number of objects in a given mass range.

    5. MassFunction           : the mass function

    6. nufnu_PS               : Finds the function nu x f(nu), used in the mass function, for Press-Schechter

    7. nufnu_ST               : Finds the function nu x f(nu), used in the mass function, for Sheth-Tormen

    8. nufnu_Jenkins          : Finds the function nu x f(nu), used in the mass function, for Jenkins

    9. nufnu_Warren           : Finds the function nu x f(nu), used in the mass function, for Warren

    10. nufnu_Reed03          : Finds the function nu x f(nu), used in the mass function, for Reed (2003)

    11. nufnu_Reed07          : Finds the function nu x f(nu), used in the mass function, for Reed (2007)

    12. N_eff                 : The effective slope of the power spectrum at the scale of the mass in question

    13. nufnu                 : The function nu x f(nu), compiled into one array for different approaches.
    """
    
    
    def __init__(self,k,Transfer,z = 0.0, step_size = 0.001 ,
                 cosmology = {"sigma_8":0.812,"n":1,"mean_dens":0.273*2.7755*10**11,"crit_dens":1.686,"hubble":0.7}, 
                 WDM = None, k_bounds = []):
        """
        Initializes the cosmology for which to perform the perturbation analysis.      
        """
        # We save these values for reference throughout the class.
        self.k = np.array(k)
        self.T = np.array(Transfer)
        self.steps = step_size
        
        self.z = z
        
        self.sigma_8 = cosmology["sigma_8"]
        self.n = cosmology["n"]
        self.mean_dens = cosmology["mean_dens"]
        self.crit_dens = cosmology["crit_dens"]
        self.hubble = cosmology["hubble"]
        
        self.WDM = WDM
            
            
        # Save the true k_bounds.
        if k_bounds:
            self.k_lower = np.log(k_bounds[0])
            self.k_high = np.log(k_bounds[1])
        else:
            self.k_lower = min(k)
            self.k_high = max(k)

        # Make sure the step-size is appropriate for romberg integration.
        self.StepSize()
        
        
        # Preliminary operations on the transfer function performed, to transform it
        # into a power spectrum that is normalized and extrapolated to limits defined
        # in config.
        self.Interpolate(k,Transfer,self.steps)
        self.Extrapolate(k,self.steps)

        if self.WDM:
            self.WDM_PS(self.WDM)
        self.Normalize(self.sigma_8)
        
 
    def StepSize(self):
        """
        Defines the step size of the resultant integration more precisely for use in Romberg integration.
        
        Input: k: the vector of ln(k) values from the transfer function.
        Output: step_size: the step-size calculated for use in the given k-bounds (in config), in Romberg Integration
        
        Note: for Romberg integration, the number of steps must be 2**p+1 where p is an integer, which is why this scaling
                should be performed.
        """
             
        index = np.ceil((np.log2((self.k_high-self.k_lower)/self.steps))-1.0)
        step_size = (self.k_high-self.k_lower)/(2**index+1)
        self.vector_length = 2**index+1
        

        
        self.steps = step_size
    def Interpolate(self,k,Transfer,step_size):
        """
        Interpolates the given Transfer function and transforms it into a Power Spectrum.
        
        Input: k        : an array of values of lnk
               Transfer : an array of values of lnT
               step_size: the step_size between values of lnk required in interpolation.
                  
        Notes: the power spectrum is calculated using P(k) =  A.k^n.T^2, with A=1 in this method
                (which leaves the power spectrum un-normalized). The interpolation is done using
                cubic splines.
        """

        # Use cubic spline interpolation over the range of values in the vector transfer function.
        interpolation = interp1d(k,config.n*k+2.0*Transfer,'cubic')
        
        # Evaluate the power spectrum (un-normalized) at new grid-spacing defined by step_size.
        self.power_spectrum = interpolation(np.arange(min(k),max(k),step_size))
         

    
    def Extrapolate(self,k,step_size):
        """
        Assumes linearity in log-log space for last 100 points and extrapolates to the extent defined in config.
        
        NOTE: There may be a more sure way of extrapolating, but this does quite alright.
        """

        # Define a new vector of wavenumbers with the new spacing.
        interp_k = np.arange(min(k),max(k),step_size)
        

        # We only extrapolate the power spectrum at the high end if the bounds are above the original bounds.
        if self.k_high > max(k):
            # Calculate slope and displacement to form a linear extrapolation
            upper_slope = (self.power_spectrum[-50]-self.power_spectrum[-1])/(interp_k[-50]-interp_k[-1])
            upper_displacement = self.power_spectrum[-1] - upper_slope*interp_k[-1]
            upper_extrapolation = upper_slope*np.arange(max(k),self.k_high,step_size)+upper_displacement
            # Join the interpolated power spectrum on to the high-k extrapolation.
            self.power_spectrum = np.concatenate((self.power_spectrum,upper_extrapolation))
        else:
            for i,ps in enumerate(self.power_spectrum):
                if ps > self.k_high:
                    high_index = i-1
                    break
               
            self.power_spectrum = self.power_spectrum[0:high_index]     

        # Only extrapolate the power spectrum at the low end if the bounds are below the original.
        if self.k_lower < min(k):
            # Calculate slope and displacement to form linear extrapolation.
            lower_slope = (self.power_spectrum[0]-self.power_spectrum[50])/(interp_k[-0]-interp_k[50])
            lower_displacement = self.power_spectrum[0] - lower_slope*interp_k[0]
            lower_extrapolation = lower_slope*np.arange(self.k_lower,min(k),step_size)+lower_displacement
            # Join the previously extrapolated (at high end if option is True) power spectrum to low-k extrapolation.
            self.power_spectrum = np.concatenate((lower_extrapolation,self.power_spectrum))
        else:
            for i,ps in enumerate(self.power_spectrum):
                if ps > self.k_lower:
                    low_index = i
                    break
               
            self.power_spectrum = self.power_spectrum[low_index:]
        
        # Form new k-vector based on new bounds and the step_size.      
        k_extended = np.arange(self.k_lower,self.k_high,step_size)   
        
        # In case there is slight overflow with vector length, restrict the vectors to the vector length
        # defined by setup.StepSize(), which needs to be precise for romberg integration.     
        self.power_spectrum = self.power_spectrum[0:self.vector_length]
        self.k_extended = k_extended[0:self.vector_length]
        

    
    def Normalize(self,norm_sigma_8):
        """
        Normalize the power spectrum to a given sigma_8
        
        Input: norm_sigma_8: The sigma_8 value to normalize to.
        """

        # Calculate the value of sigma_8 without prior normalization.

        sigma_8 = self.MassVariance((4.*np.pi*8**3*self.mean_dens/3.))
        
        # Calculate the normalization factor A.
        self.normalization = norm_sigma_8**2/sigma_8
    
        # Normalize the previously calculated power spectrum.
        self.power_spectrum = np.log(self.normalization)+self.power_spectrum
     
    def Radius(self,M):
        """
        Calculates radius from mass given mean density.
        """
        
        self.R = (3.*M/(4.*np.pi*self.mean_dens))**(1./3.)
        
        return self.R
    def TopHat_WindowFunction(self,M):
        """
        Constructs the window function squared in Fourier space for given radii
        
        Input: R: The radius of the top-hat function
        Output: W_squared: The square of the top-hat window function in Fourier space.
        """
        # Calculate the factor kR, minding to un-log k before use.
        kR = np.exp(self.k_extended)*self.Radius(M)
        W_squared = (3*(np.sin(kR) - kR*np.cos(kR))/kR**3)**2
        

        return W_squared
    

    def MassVariance(self,M):
        """
        Finds the Mass Variance of R using the top-hat window function.
        
        Input: R: the radius of the top-hat function (single number).
        Output: sigma_squared: the mass variance.
        """
        # Calculate the integrand of the function. Note that the power spectrum and k values must be
        #'un-logged' before use, and we multiply by k because our steps are in logk. 
        integrand = np.exp(self.k_extended)**3*np.exp(self.power_spectrum)*self.TopHat_WindowFunction(M)
                
        # Perform the integration using romberg integration, and multiply be the pre-factors.
        sigma_squared = (0.5/np.pi**2)*intg.romb(integrand,dx = self.steps)
                
                
        # Multiply by the growth factor for the redshifts involved.
        dist = Distances(camb_config.omega_lambda, camb_config.omega_baryon+camb_config.omega_cdm,0.0,H_0 = camb_config.hubble)
        growth = dist.GrowthFactor(self.z)**2
       
        sigma_squared = sigma_squared*growth
        
        return sigma_squared
    
     
    def dlnsdlnM(self,M):
        """
        Uses a top-hat window function to calculate |dlnsigma/dlnM| at a given radius.
        
        Input: R: the radius of the top-hat function. 
        Output: integral: the derivatiave of log sigma with respect to log Mass.
        """
        

        R = self.Radius(M)
        # define the vector g = kR
        g = np.exp(self.k_extended)*R
            
        # Find the mass variance.
        sigma = np.sqrt(self.MassVariance(M))
            
        # Define the derivative of the window function squared.
        window = (np.sin(g)-g*np.cos(g))*(np.sin(g)*(1-3.0/(g**2))+3.0*np.cos(g)/g)
            
        # Compile into the integrand.
        integrand = (np.exp(self.power_spectrum)/np.exp(self.k_extended))*window
            
        # Integrate using romberg integration and multiply by pre-factors.
        integral =(3.0/(2.0*sigma**2*np.pi**2*R**4))*intg.romb(integrand,dx=self.steps)
        
        return integral
    
    def dndM(self,M,approach = "PS"):
        """
        Computes the value of dn/dM for a given radius.
        
        Input: R: radius of the top-hat function.
              Approach: A string indicating which approach to take to the fitting function. 
                        Valid values are:
                        1. PS: Press-Schechter Approach
                        2. ST: Sheth-Tormen
                        3. Jenkins: Jenkins empirical fit
                        4. Warren: Warren empirical fit
                        5. Reed03: Reed empirical from 2003
                        6. Reed07: Reed empirical from 2007
                        
                      All fits are taken from Lukic et. al. 2007
        Output dndm: the number density of objects in a mass range dM.
        """

        
        leftovers = config.mean_dens/M**2
        
        dndm = self.nufnu(M,approach) *leftovers*np.abs(self.dlnsdlnM(M))
        
        return dndm
    
    def MassFunction(self,M,approach = "PS"):
        """
        Uses the Press Schechter approach with a spherical top-hat to calculate the mass function.
        
        Input: R: radius of the top-hat function
              Approach: A string indicating which approach to take to the fitting function. 
                        Valid values are:
                        1. PS: Press-Schechter Approach
                        2. ST: Sheth-Tormen
                        3. Jenkins: Jenkins empirical fit
                        4. Warren: Warren empirical fit
                        5. Reed03: Reed empirical from 2003
                        6. Reed07: Reed empirical from 2007
                        
                      All fits are taken from Lukic et. al. 2007
        Output mass_function: the Press-Schechter mass function log10(dn/dlnM).
        """

        mass_function = np.log10(np.log(10.0)*M*self.dndM(M,approach))

        return mass_function
    
    def nufnu_PS(self,M):
        """
        Computes the function nu.f(nu) for the Press-Schechter approach at a given radius.
        
        Input R: radius of the top-hat function
        Output: f_of_nu: the function nu.f(nu) for the PS approach.
        """
        

        sigma = np.sqrt(self.MassVariance(M))
        f_of_nu = np.sqrt(2.0/np.pi)*(self.crit_dens/sigma)*np.exp(-0.5*(self.crit_dens/sigma)**2)   
        
        return f_of_nu    
    
    def nufnu_ST(self,M):
        """
        Finds the Sheth Tormen vf(v) 
        
        Input R: radius of the top-hat function
        Output: vfv: the Sheth-Tormen mass function fit.
        """

        nu = config.crit_dens/np.sqrt(self.MassVariance(M))
        a = 0.707
        
        vfv = 0.3222*np.sqrt(2.0*a/np.pi)*nu*np.exp(-(a*nu**2)/2.0)*(1+(1.0/(a*nu**2))**0.3)
        
        return vfv
    
    def nufnu_Jenkins(self,M):
        """
        Finds the Jenkins empirical vf(v) 
        
        Input R: radius of the top-hat function
        Output: vfv: the Jenkins mass function fit.
        """

        sigma = np.sqrt(self.MassVariance(M))
        
        vfv = 0.315*np.exp(-np.abs(np.log(1.0/sigma)+0.61)**3.8)
        
        return vfv
    
    def nufnu_Warren(self,M):
        """
        Finds the Warren empirical vf(v) 
        
        Input R: radius of the top-hat function
        Output: vfv: the Warren mass function fit.
        """

        sigma = np.sqrt(self.MassVariance(M))
        
        vfv = 0.7234*((1.0/sigma)**1.625+0.2538)*np.exp(-1.1982/sigma**2)
        
        return vfv
    
    def nufnu_Reed03(self,M):
        """
        Finds the Reed 2003 empirical vf(v) 
        
        Input R: radius of the top-hat function
        Output: vfv: the Reed 2003 mass function fit.
        
        NOTE: Only valid from -1.7 < ln sigma^-1 < 0.9
        """
        
   
        sigma = np.sqrt(self.MassVariance(M))
                
        ST_Fit = self.nufnu_ST(M)
        
        vfv = ST_Fit*np.exp(-0.7/(sigma*np.cosh(2.0*sigma)**5))
        
        return vfv
    
    def nufnu_Reed07(self,M):
        """
        Finds the Reed 2003 empirical vf(v) 
        
        Input R: radius of the top-hat function
        Output: vfv: the Reed 2003 mass function fit.
        
        NOTE: Only valid from -1.7 < ln sigma^-1 < 0.9
        """
        sigma = np.sqrt(self.MassVariance(M))
        nu = config.crit_dens/sigma
        
        G_1 = ((0.4-1.0/sigma)**2)**(0.5/0.6**2)
        G_2 = ((0.75-1.0/sigma)**2)**(0.5/0.2**2)
        
        c = 1.08
        a = 0.764/c
        A = 0.3222
        p = 0.3
        n_eff = self.N_Eff(M)
        
        if sigma > np.exp(-0.55) and sigma < np.exp(1.7):
            vfv = A*np.sqrt(2.0*a/np.pi)*(1.0+(1.0/(a*nu**2))**p+0.6*G_1+0.4*G_2)*nu*np.exp(-c*a*nu**2/2.0-0.03*nu**0.6/(n_eff+3)**2)
        else:
            vfv = np.NaN
        return vfv

    def N_Eff(self,M):
        """
        Calculates the power spectral slope at the scale of the halo radius, using eq. 42 in Lukic et. al 2007.
        
        """
        
        n_eff = -3.0 *(2.0*self.dlnsdlnM(M)+1.0)
        
        return n_eff
        
    def WDM_PS(self,WDM):
        """
        Tansforms the CDM Power Spectrum into a WDM power spectrum for a given warm particle mass m_x.
        
        NOTE: formula from Bode et. al. 2001 eq. A9
        """

        h = self.hubble
        g_x = 1.5
        m_x = WDM
        nu = 1.12
            
        alpha = 0.049*(camb_config.omega_cdm/0.25)**0.11*(h/0.7)**1.22*(1/m_x)**1.11*(1.5/g_x)**0.29
            
        Transfer = (1+(alpha*np.exp(self.k_extended))**(2*nu))**-(5.0/nu)
        trans = open(config.out_pref+"/WDM"+str(WDM)+"Transfer.txt",'w')
        for i,t in enumerate(Transfer):
            trans.write(str(self.k_extended[i])+"     "+str(t)+"\n")
                
        power_spectrum = np.log(np.exp(self.power_spectrum)*Transfer**2)

        self.power_spectrum = power_spectrum
        
    def nufnu(self,M,approach="PS"):
        """
        Sets the nufnu based on the approach
        """
        if approach is "PS":
            nufnu = self.nufnu_PS(M) 
        elif approach is "ST":
            nufnu = self.nufnu_ST(M)
        elif approach is "Jenkins":
            nufnu = self.nufnu_Jenkins(M)
        elif approach is "Warren":
            nufnu = self.nufnu_Warren(M)
        elif approach is "Reed03":
            nufnu = self.nufnu_Reed03(M) 
        elif approach is "Reed07":
            nufnu = self.nufnu_Reed07(M) 
        else:
            sys.exit("Approach has invalid argument")
        
        return nufnu

class Distances(object):
    """
    A Class to calculate distances and other cosmological quantities.
    """        
    
    def __init__(self,omega_lambda = 0.7,omega_mass = 0.3,omega_k = 0.0,H_0 = 70):
        
        self.omega_lambda = omega_lambda
        self.omega_mass = omega_mass
        self.omega_k = omega_k
        self.H_0 = H_0
        
    def ScaleFactor(self,z):
        """
        Calculates the scale factor at a given redshift
        
        Input: z: the redshift.
        Output a: the scale factor at the given redshift.
        """
        
        a = 1.0/(1.0+z)

        return a
    
    def zofa(self,a):
        """
        Returns the redshift of a particular scale factor
        """
        
        z = 1.0/a -1.0
        
        return z
    
    def HubbleFunction(self,z):
        """
        Finds the hubble parameter at z normalized by the current hubble constant.
        
        Input: z: redshift
        Output: hubble: the hubble parameter at z divided by H_0
        """
        
        hubble = np.sqrt(self.omega_mass/self.ScaleFactor(z)**3 + (1.0-self.omega_mass))
        
        return hubble
    
    def Dplus(self,redshifts):
        """
        Finds the factor D+(a), from Lukic et. al. 2007, eq. 8.
        
        Uses romberg integration with a suitable step-size. 
        
        Input: z: redshift.
        
        Output: dplus: the factor.
        """
        if type(redshifts) is type (0.4):
            redshifts = [redshifts]
            
        dplus = np.zeros_like(redshifts)
        for i,z in enumerate(redshifts):
            step_size = self.StepSize(0.0000001,self.ScaleFactor(z))
            a_vector = np.arange(0.0000001,self.ScaleFactor(z),step_size)
            integrand = 1.0/(a_vector*self.HubbleFunction(self.zofa(a_vector)))**3
            
            integral = intg.romb(integrand,dx=step_size)
            dplus[i] = 5.0*self.omega_mass*self.HubbleFunction(z)*integral/2.0
        
        return dplus
    
    def GrowthFactor(self,z):
        """
        Finds the factor d(a) = D+(a)/D+(a=1), from Lukic et. al. 2007, eq. 8.
        
        Input: z: redshift.
        
        Output: growth: the growth factor.
        """
        
        growth = self.Dplus(z)/self.Dplus(0.0)

        return growth
    
    def StepSize(self,mini,maxi):
        """
        Calculates a suitable step size for romberg integration given data limits
        """
        
        p = 13
        
        while (maxi-mini)/(2**p+1) < 10**(-5):
            p=p-1
            
        step_size = (maxi-mini)/(2**p+1.0)
        
        return step_size
        
        

        
        
    