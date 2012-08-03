'''
Created on Jun 12, 2012

@author: Steven
'''
from config import *

###################################################################################################
if __name__ == '__main__':
    
    from PyCosmology.sims.load import sims
    
    sim = sims(sim_file,group_filename,ids_filename,n_groups,reverse,centring_mode,loud=loud,subgroup_filename=subgroup_file)
    
    if plot_density_profile:
        sim.plot_group_density_profile()
    if plot_group_pos:
        sim.plot_group_pos()
    if plot_vel_disp:
        sim.plot_velocity_disp()
    if plot_shell_structure:
        sim.structure_in_shells(radial_bins)
    
    if plot_axis_ratios:
        sim.AxisRatioDistribution()
    if plot_mass_offset:
        sim.Centre_of_mass_offsets()
    if plot_non_param:
        sim.non_parametric_shell_structure(n_radial_bins,n_angles,nside)