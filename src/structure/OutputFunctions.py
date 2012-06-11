'''
Created on Apr 20, 2012

@author: jeven
'''
import matplotlib
matplotlib.use("AGG")
from matplotlib import pyplot as plt 
import numpy as np
import TableIO


from Config import *

def PowerPlot(k,ps):
    """
    Output a png image of the plot of the power spectrum after interpolation and extrapolation
    """
    print "  Producing Power Spectrum Plot"
    
 
    plt.clf()
    plt.title('Logarithmic Matter Power Spectrum')
    plt.grid(True)  
    plt.xlabel(r'Wavenumber, $\log(k)$')  
    plt.ylabel(r'Power, $\log[P(k)]$') 
    
    linetypes = ['r-','b-','g-','m-','k-']
    if len(config.WDM)+1 > len(linetypes):
        print "WARNING: more power spectra (", len(config.WDM),") than linetypes (", len(linetypes),"), some will be used twice."
    for i,wdm in enumerate(ps):
        if i ==0:
            plt.plot(k,wdm, linetypes[i], lw=1, label = "CDM Model")
        else:
            plt.plot(k,wdm, linetypes[i], lw=1, label = "WDM Particle Mass "+str(config.WDM[i-1])+"keV")
    plt.legend(loc=0)
                                
    plt.savefig(config.out_pref+'/PowerSpectrum.eps')
        
def SigmaPlot(M,sigma):
    """
    Output a png image of the plot of sigma for a range of R
    """
    print "  Producing Mass Variance Plot"
        
    plt.clf()
    plt.title('Mass Variance')
    plt.grid(True)  
    plt.xlabel(r'Mass $(\log_{10} M_{sun})$')  
    plt.ylabel(r'$\sigma (M)$') 
    
    linecolours = ['r','b','g','m','k'] # associated with WDM sizes
    linetypes = ['-','o','--','v','s','*','^'] # associated with redshifts
    
    if len(config.WDM)+1 > len(linecolours):
        print "WARNING: more wdm particles (", len(config.WDM),") than line colours (", len(linecolours),"), some will be used twice."
    if len(config.z) > len(linetypes):
        print "WARNING: more redshifts (", len(config.WDM),") than line types (", len(linetypes),"), some will be used twice."
    
    for i in range(sigma.shape[0]):
        for j in range(sigma.shape[1]):
            if j == 0:
                plt.plot(M,sigma[i,j,:], linecolours[j]+linetypes[i],lw=1,label = "CDM Model, Redshift = "+str(config.z[i])  )
            else:
                plt.plot(M,sigma[i,j,:], linecolours[j]+linetypes[i],lw=1,label = "WDM Particle Mass "+str(config.WDM[j-1])+"keV, Redshift = "+str(config.z[i])  )
    plt.legend(loc=0)                            
    plt.savefig(config.out_pref+'/MassVariance.eps')


def SigmaLogPlot(M, sigma):
    """
    Output a png image of the plot of log sigma squared for a range of R
    """
    print "  Producing Logarithmic Mass Variance Plot"
        
    plt.clf()
    plt.title('Logarithmic Mass Variance')
    plt.grid(True)  
    plt.xlabel(r'Mass $(\log_{10} M_{sun})$')  
    plt.ylabel(r'$\log[ \sigma^2 (R)]$') 
    linecolours = ['r','b','g','m','k'] # associated with WDM sizes
    linetypes = ['-','o','--','v','s','*','^'] # associated with redshifts
    
    if len(config.WDM)+1 > len(linecolours):
        print "WARNING: more wdm particles (", len(config.WDM),") than line colours (", len(linecolours),"), some will be used twice."
    if len(config.z) > len(linetypes):
        print "WARNING: more redshifts (", len(config.z),") than line types (", len(linetypes),"), some will be used twice."
    
    for i in range(sigma.shape[0]):
        for j in range(sigma.shape[1]):
            if j == 0:
                plt.plot(M,np.log10(sigma[i,j,:]), linecolours[j]+linetypes[i],lw=1,label = "CDM Model, Redshift = "+str(config.z[i])  )
            else:
                plt.plot(M,np.log10(sigma[i,j,:]), linecolours[j]+linetypes[i],lw=1,label = "WDM Particle Mass "+str(config.WDM[j-1])+"keV, Redshift = "+str(config.z[i])  )
    plt.legend(loc=0)                           
    plt.savefig(config.out_pref+'/LogMassVariance.eps')


def StepSizePlot(steps,sigma,prefix):
    """
    Plots sigma_8 for various step sizes of the integration.
    """
    print "  Producing Plot of Step Size Effect"
        
    sigma_8_true = []
    for i in range(0,config.step_size_number):
        sigma_8_true = sigma_8_true+[config.sigma_8]

    plt.clf()
    plt.title(r'$\sigma_8$ for various step-sizes of integration.')
    plt.grid(True)  
    plt.xlabel(r'Step Size $\log_{10}$')  
    plt.ylabel(r'$\sigma_8$') 
    plt.plot(np.log10(steps),np.sqrt(sigma), 'ro')
    plt.plot(np.log10(steps),sigma_8_true,'b-')
                                
    plt.savefig(config.out_pref+'/StepSizeEffect.png')

      
def MassFunctionPlot(mass,mf,sims=None,mass_bins=None):
    """
    Writes an image of the mass function
    """
    print "  Producing Plot of Mass Function"
        
    plt.clf()
    plt.figure(figsize=(10,10))
    plt.title('Mass Function')
    plt.grid(True)  
    plt.xlabel(r'Mass $(\log M_{sun})$')  
    plt.ylabel(r'Logarithmic Mass Function $\log \left( \frac{dn}{d \ln M} \right) $') 
    
    linecolours = ['r','b','g','m','k','c'] # associated with approaches
    linetypes = ['-','-.',':','.','--','v','s','*','^','o'] # associated with WDM sizes and redshifts

    
    if len(config.WDM)+1 > len(linecolours):
        print "WARNING: more wdm particles (", len(config.WDM),") than line colours (", len(linecolours),"), some will be used twice."
    if len(config.z) > len(linetypes):
        print "WARNING: more redshifts (", len(config.z),") than line types (", len(linetypes),"), some will be used twice."
    if len(config.approach) >len(linetypes):
        print "WARNING: more redshifts (", len(config.approach),") than line types (", len(linetypes),"), some will be used twice."
    
    for i in range(mf.shape[0]):
        for j in range(mf.shape[1]):
            for k in range(mf.shape[2]):
                        
                if k == 0:
                    plt.plot(mass,mf[i,j,k,:], linecolours[i]+linetypes[mf.shape[2]*j+k],lw=1,label = "CDM, "+config.approach[i]+", z = "+str(config.z[j])  )
                else:
                    plt.plot(mass,mf[i,j,k,:], linecolours[i]+linetypes[mf.shape[2]*j+k],lw=1,label = "WDM "+str(config.WDM[k-1])+"keV,"+config.approach[i]+" z = "+str(config.z[j])  )
    
          
    if sims:
        for i,sim in enumerate(sims):
            plt.plot(mass_bins[i],sim, linecolours[i+mf.shape[0]+1]+linetypes[0],lw=1,label=config.sim_files[i])
            
    plt.legend(loc=0)                            
    plt.savefig(config.out_pref+'/MassFunction.eps')
        
        
def EverythingOutput(data, z, approach, wdm="CDM"):
    """
    Outputs most of the calculated values as a .dat file
    
    Input: data: a list of lists. An array containing all of the data to write to a table.
    """
    print "  Outputting Data Tables"
    if wdm is None:
        wdm = "CDM"
    mass_var_table = TableIO.TableFromColumns(data)
    if wdm is "CDM":
        mass_var_table['filename'] = config.out_pref+'/OutputTables/'+'z'+str(z)+'.CDM.'+approach+'.txt'
    else:
        mass_var_table['filename'] = config.out_pref+'/OutputTables/'+'z'+str(z)+'.WDM'+str(wdm)+'.'+approach+'.txt'
    TableIO.writeTable(mass_var_table)