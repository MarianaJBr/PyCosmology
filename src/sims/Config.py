'''
Created on May 8, 2012

@author: smurray
'''

class Config(object):
    '''
    A simple config file for the main simulaiton processing script.
    '''

    sim_directory = "/Users/smurray/Documents/PhD/Simulations/data_for_steven/"
    
    pos_file = "s130_pos_zslice.dat"
    group_pos_file = "groups/GroupData/Group_1/s130_pos_rotated.dat"
    group_density_file = "groups/GroupData/Group_1/s130_densities.dat"
    group_velocity_file = "groups/GroupData/Group_1/s130_vel.dat"
    group_velocity_r_file = "groups/GroupData/Group_1/s130_vel_rotated.dat"
    dot_plot_file = "density_plot.eps"
    density_file = "DensityProfile.eps"
    vd_original_file = "groups/GroupData/Group_1/VelocityDispersionOriginal.eps"
    vd_rotated_file = "groups/GroupData/Group_1/VelocityDispersionRotated.eps"
    
    
    density_smoothing = "sgolay"
    