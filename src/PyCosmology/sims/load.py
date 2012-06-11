'''
Created on May 7, 2012

@author: smurray
'''
import numpy as np

from sims.fort.ReadSim import simops as fsims
import plotting

import os, errno

class sims(object):
    """
    A class to import simulation data and group data from gadget runs. Basic plotting/manipulation.
    """
    
    def __init__(self,filename,group_filename=None,ids_filename=None,n_groups=1,reverse=False,centring_mode=None):
        """
        Imports the simulation and detects groups if specified.
        """
    
        # readsim reads in the simulation and defines fsims.pos and fsims.vel amongst others.
        fsims.readsim(filename)

        # find_groups locates groups in the simulation and saves the array fsmims.group_number,
        # which contains the number of the group for each particle.  
        if group_filename:
            fsims.find_groups(group_filename,ids_filename)
            
            self.groups_pos = []
            self.groups_vel = []
            for i in range(n_groups):
                if not reverse:
                    self.groups_pos = self.groups_pos + [np.asfortranarray(fsims.pos[:,fsims.group_number == i+1])]
                    self.groups_vel = self.groups_vel + [np.asfortranarray(fsims.vel[:,fsims.group_number == i+1])]
                else:
                    self.groups_pos = self.groups_pos + [np.asfortranarray(fsims.pos[:,fsims.group_number == np.max(fsims.group_number)-i])]
                    self.groups_vel = self.groups_vel + [np.asfortranarray(fsims.vel[:,fsims.group_number == np.max(fsims.group_number)-i])]
            
                fsims.centre(self.groups_pos[i],mode=centring_mode)
                fsims.rotate(self.groups_pos[i],self.groups_vel[i],[0.0,0.0,0.0],0.0)    

        self.make_directories(filename)
        
    def make_directories(self,filename):
        """
        Makes a simple directory structure to store images etc.
        """
        self.sim_dir = filename.partition('.')[0]
        try:
            os.makedirs(self.sim_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
        print self.sim_dir    
        os.chdir(self.sim_dir)
        
        self.group_dir = "GroupData"
        try:
            os.makedirs(self.group_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
        
        self.group_density_dir = "GroupData/DensityPlots"
        try:
            os.makedirs(self.group_density_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
            
        self.group_profile_dir = "GroupData/DensityProfiles"
        try:
            os.makedirs(self.group_profile_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
        self.group_velocity_dir = "GroupData/VelocityDispersion"
        try:
            os.makedirs(self.group_velocity_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise        
            
    def plot_group_pos(self):
        """
        Convenience function for plotting the group positions as density maps all at once.
        """
        
        for i,group in enumerate(self.groups_pos):
            plotting.SpatialPlot(group[0,:],group[1,:],self.group_density_dir+'/group_'+str(i),smoothing='sgolay')
            
    def plot_group_density_profile(self):
        """
        Convenience function for plotting the density profiles for all maintained groups.
        """
        
        for i,group in enumerate(self.groups_pos):
            dens,edges = fsims.density_profile(group,10,0.01)
            plotting.DensityProfile(edges[:-1],dens,self.group_profile_dir+'/group_'+str(i))
        
    def plot_velocity_disp(self):
        """
        Convenience function for plotting the velocity dispersion profiles for all maintained groups.
        """
        
        for i,group in enumerate(self.groups_vel):
            plotting.VelocityDispersion(group.T, self.group_velocity_dir+'/group_'+str(i))

