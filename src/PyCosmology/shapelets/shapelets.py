'''
Created on Jun 5, 2012

@author: Steven
'''

from PyCosmology.sims.fort.read_sim import simops
from PyCosmology.shapelets.fort.shapelets import shapelets as fc
import numpy as np
from scipy.misc import comb
#import matplotlib
#matplotlib.use("Qt4Agg")
#import matplotlib.pyplot as plt


class ShapeletOperations_DMH(object):
    """
    This class contains methods for evaluating shapelet co-efficients and plotting
    results based on this.
    """
    
    def __init__(self,sim_file,groups_file,group_ids_file, n=1,reverse=False):
        """
        Read in a sim-file for use with the shapelet analysis.
        
        PARAMETERS:
        sim_file: file containing the simulation snapshot
        groups_file: file containing the location of groups
        group_ids_file: file containing the ids of the group members
        
        n: optional, default 1, number of groups to retain for further analysis (from largest to smallest)
        reverse: optional, default False, whether to use the n SMALLEST groups rather than largest.
        """
        
        simops.readsim(sim_file)
        simops.find_groups(groups_file,group_ids_file)
        
        self.groups_pos = []
        for i in range(n):
            if not reverse:
                self.groups_pos = self.groups_pos + [np.asfortranarray(simops.pos[:,simops.group_number == i+1])]
            else:
                self.groups_pos = self.groups_pos + [np.asfortranarray(simops.pos[:,simops.group_number == np.max(simops.group_number)-i])]
            
            simops.centre(self.groups_pos[i],mode=3)
            simops.rotate(self.groups_pos[i],self.groups_pos[i],[0.0,0.0,0.0],0.0)

        
    def find_coeffs(self,bins):
        """
        Finds the coefficients for the shapelets
        """
        self.coefficients = []
        self.beta = []
        self.nmax = (bins-3)/2
        
        for i,group in enumerate(self.groups_pos):
            x_max = 2.0*np.max(group[0,:])

            hist,edges = np.histogramdd(group.T,bins=bins,range=[[-x_max,x_max],[-x_max,x_max],[-x_max,x_max]])
            del edges
            print np.sum(hist)

            hist = simops.massarr[1]*hist

            self.coefficients = self.coefficients + [np.asfortranarray(fc.coeff_cube(np.asfortranarray(hist),x_max))]
            self.beta = self.beta+[np.max(group[0,:])/np.sqrt(2.*bins)]
            
            print np.sum(self.coefficients[i])
            print self.coefficients[i]
    def reconstruct(self,bins):
        """
        Reconstructs the 'original' data using the basis functions and amplitudes
        """
        
        self.rebuilt = []
        for i,coefficients in enumerate(self.coefficients):
            self.rebuilt = self.rebuilt + [fc.reconstruct_cube(bins,coefficients,(2.0*np.min(self.groups_pos[i][0,:]),2.0*np.max(self.groups_pos[i][0,:])))]
            
                             
    def zeroth_moment(self):
        """
        Finds the zeroth moment of the shape based on the coefficients found. Associated with total mass.
        """
        self.M_0 = []
        a = np.arange(0,self.nmax,2)
        print a
        
        for i,coeff in enumerate(self.coefficients):
            self.M_0 = self.M_0+[0.0]

            for n1 in a:
                for n2 in a:
                    for n3 in a:
                        self.M_0[i] = self.M_0[i]  + coeff[n1,n2,n3]*self.U(n1+n2+n3)*self.W(n1,n2,n3)
            
            
            self.M_0[i] = np.sqrt(np.sqrt(np.pi))**3 *np.sqrt(self.beta[i])**3
            
        print "zeroth order moment (total mass): ", self.M_0[i]
        
    def U(self,n):
        """
        Calculates an instance of U from Fluke et. al
        """
        return 2**((3-n)/2)
    
    def W(self,n1,n2,n3):
        """
        Calculates an instance of W from Fluke. et. al.
        """
        
        return np.sqrt(comb(n1,n1/2)*comb(n2,n2/2)*comb(n3,n3/2))
        
    def J(self,n,beta):
        """
        Calculates an instance of J from Fluke et. al.
        """
        
        return np.sqrt(2**(1-n)*np.sqrt(np.pi)*beta)*np.sqrt(comb(n,n/2))
        
    def Centroid(self):
        """
        Finds the centroid of a group defined by Fluke et. al.
        
        Hmm..
        """
        
        self.centroid = []
        a = np.arange(0,self.nmax,2)
        b = np.arange(1,self.nmax,2)
        
        for i,coeff in enumerate(self.coefficients):
            x1,y1,z1 = (0.0,0.0,0.0)
            for n1 in b:
                for n2 in a:
                    for n3 in a:
                        x1 = x1 +coeff[n1,n2,n3]*np.sqrt(n1+1)*self.U(n1+n2+n3)*self.W(n1+1,n2,n3)
                        
            x1 = x1*np.pi**(3./4.)*self.beta[i]**(5./2.)/self.M_0[i]
             
            for n2 in b:
                for n1 in a:
                    for n3 in a:
                        y1 = y1+coeff[n1,n2,n3]*np.sqrt(n2+1)*self.U(n1+n2+n3)*self.W(n1,n2+1,n3)
                
            y1 = y1*np.pi**(3./4.)*self.beta[i]**(5./2.)/self.M_0[i]
            
            for n3 in b:
                for n1 in a:
                    for n2 in a:
                        z1 = z1+coeff[n1,n2,n3]*np.sqrt(n3+1)*self.U(n1+n2+n3)*self.W(n1,n2,n3+1)
                
            z1 = z1*np.pi**(3./4.)*self.beta[i]**(5./2.)/self.M_0[i]
        
            self.centroid = self.centroid = [(x1,y1,z1)]
            print self.centroid[i]
                    
            
        