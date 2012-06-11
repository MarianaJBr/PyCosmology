'''
Created on May 7, 2012

@author: smurray
'''
import TableIO
import numpy as np
from scipy.optimize import curve_fit

import PlottingRoutines as pr

from Config import *

if __name__ == '__main__':
    
    x,y,z = TableIO.readColumns(Config.sim_directory+Config.pos_file, "!#",[0,1,2])
    
    pr.SpatialPlot(x,y,Config.sim_directory+Config.dot_plot_file,smoothing=Config.density_smoothing,
                   condition = 'with 9.8 < z < 10.2',smoothing_scale=0.05,resolution = 400)
    
    
    
    group_x, group_y, group_z = TableIO.readColumns(Config.sim_directory+Config.group_pos_file, "!#",[0,1,2])
    pr.SpatialPlot(group_x,group_y,Config.sim_directory+'groups/GroupData/Group_1/'+Config.dot_plot_file,smoothing=Config.density_smoothing,
                   condition = 'for the largest halo',smoothing_scale=0.05,resolution = 400)
    
    
    
    r,dens = TableIO.readColumns(Config.sim_directory+Config.group_density_file,"!#",[0,1])
    
    
    def NFW(r,rho_0, Rs):
        """
        The NFW profile in terms of its variables
        """
        return np.log(rho_0/((r/Rs)*(1.0+r/Rs)**2))
    
    
    #popt, pcov = curve_fit(NFW, r, dens,p0=(5000.0,0.5))
    
    #print "Optimal Value of rho_0: ", popt[0]
    #print "Variance of rho_0: ", pcov[0,0]
    
    #print "Optimal Value of Rs: ", popt[1]
    #print "Variance of Rs: ", pcov[1,1]
    
    #r_NFW = np.linspace(min(r), max(r), 100)
    #dens_NFW = np.zeros_like(r_NFW)
    #for i,r_n in enumerate(r_NFW):
    #    dens_NFW[i] = NFW(r_n,5000,0.7)
        
    pr.DensityProfile(r,dens,Config.sim_directory+'groups/GroupData/Group_1/'+Config.density_file)
    
    
    vel = TableIO.readColumns(Config.sim_directory+Config.group_velocity_file,"!#",[0,1,2])
    pr.VelocityDispersion(vel,Config.sim_directory+Config.vd_original_file)
    vel = TableIO.readColumns(Config.sim_directory+Config.group_velocity_r_file,"!#",[0,1,2])
    pr.VelocityDispersion(vel,Config.sim_directory+Config.vd_rotated_file)
