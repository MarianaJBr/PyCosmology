'''
Created on Jun 6, 2012

@author: Steven
'''

if __name__ == '__main__':
    import shapelets
    from config import config
    import time
    
    initial = time.time()
    sh = shapelets.ShapeletOperations_DMH(config.sim_file,config.groups_file,config.ids_file,config.n_groups,config.reverse)
    
    sub_initial = time.time()
    sh.find_coeffs(config.bins)
    sub_final = time.time()
    print "  Coefficient finding time: ", sub_final - sub_initial
    
    #sub_initial = time.time()
    #sh.reconstruct(config.bins)
    #sub_final = time.time()
    #print "  Reconstruction time: ", sub_final - sub_initial
    
    sub_initial = time.time()
    sh.zeroth_moment()
    sub_final = time.time()
    print "  Zeroth Moment time: ", sub_final - sub_initial
    
    sub_initial = time.time()
    sh.Centroid()
    sub_final = time.time()
    print "  Centroid time: ", sub_final - sub_initial
    
    final = time.time()
    print "Overall Time Taken: ", final - initial