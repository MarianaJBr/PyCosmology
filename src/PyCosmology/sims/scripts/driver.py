'''
Created on Jun 12, 2012

@author: Steven
'''
# Configuration Options for the Driving Program.##################################################
sim_file = '/Users/Steven/Documents/PhD/Simulations/Snap_130/snapshot_130.cat'
group_filename = '/Users/Steven/Documents/PhD/Simulations/Snap_130/groups/groups_130.fofcat.for'
ids_filename = '/Users/Steven/Documents/PhD/Simulations/Snap_130/groups/groups_130.pos.ids'

n_groups = 100
reverse = False
centring_mode = 1

loud = False

###################################################################################################
if __name__ == '__main__':
    
    from PyCosmology.sims.load import sims
    
    sim = sims(sim_file,group_filename,ids_filename,n_groups,reverse,centring_mode,loud=loud)
    
    
    sim.plot_group_density_profile()
    sim.plot_group_pos()
    sim.plot_velocity_disp()
    
    sim.AxisRatioDistribution()