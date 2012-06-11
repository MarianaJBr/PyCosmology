'''
Created on Jun 6, 2012

@author: Steven
'''

class config(object):
    '''
    Configuration file for main.py in shapelets package.
    '''

    sim_file = "/Users/Steven/Documents/PhD/Simulations/Snap_130/snapshot_130.cat"
    groups_file = "/Users/Steven/Documents/PhD/Simulations/Snap_130/groups/groups_130.fofcat.for"
    ids_file = "/Users/Steven/Documents/PhD/Simulations/Snap_130/groups/groups_130.pos.ids"
    
    reverse = False      # If True, use smallest groups in the simulation
    n_groups = 1        # Number of groups to use.
    
    bins = 50
    