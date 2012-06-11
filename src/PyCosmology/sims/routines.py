'''
Created on Jun 11, 2012

@author: Steven
'''
import numpy as np

def NFW(r,rho_0, Rs):
    """
    The NFW profile in terms of its variables
    """
    return np.log(rho_0/((r/Rs)*(1.0+r/Rs)**2))