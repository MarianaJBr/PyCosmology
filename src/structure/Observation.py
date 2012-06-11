'''
Created on May 1, 2012

@author: smurray
'''
import numpy as np
import TableIO

class HMFObservations(object):
    """
    Performs operations on a (set of) file(s) to calculate the 'observed' HMF.
    """
    
    
    def __init__(self,filename, size):
        """
        Reads in the file and important information
        
        Input:
        file: the file with the simulation/survey results
        size: the length of a side of the simulation box
        """
        
        self.size = size
        
        
        self.masses = TableIO.readColumns(filename,"!#",columns=[1])[0]
        self.masses = np.log10(np.abs(self.masses))+10

    def MakeHist(self,bin_size):
        """
        Bins the virial mass data and divides by volume and bin width.
        """
        
        bins = np.ceil((max(self.masses)-min(self.masses))/bin_size)
        bin_size = (max(self.masses)-min(self.masses))/bins

        hist , mass_bins = np.histogram(self.masses, int(bins))

        hist = np.log10(hist/(bin_size*self.size**3))
        
        return (hist, mass_bins)
        