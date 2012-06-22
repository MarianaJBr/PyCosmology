'''
Created on Jun 6, 2012

@author: Steven
'''

if __name__ == '__main__':
    from PyCosmology.shapelets.shapelets import ShapeletOperations_DMH as shops
    from config import config
    import time
    
    initial = time.time()
    sh = shops(config.sim_file,config.groups_file,config.ids_file,config.n_groups,config.reverse,config.test,config.basis_test,config.n_test)
    
    if not config.basis_test:
        sub_initial = time.time()
        sh.find_coeffs(config.bins)
        sub_final = time.time()
        print "  Coefficient finding time: ", sub_final - sub_initial
        
        sub_initial = time.time()
        sh.zeroth_moment()
        sub_final = time.time()
        print "  Zeroth Moment time: ", sub_final - sub_initial
        
        sub_initial = time.time()
        sh.Centroid()
        sub_final = time.time()
        print "  Centroid time: ", sub_final - sub_initial
    
        sub_initial = time.time()
        sh.rmsRadius()
        sub_final = time.time()
        print "  rms Radius time: ", sub_final - sub_initial
    
        sub_initial = time.time()
        sh.InertiaTensor()
        sub_final = time.time()
        print "  Moment of Inertia time: ", sub_final - sub_initial
 
        sub_initial = time.time()
        sh.Triaxiality()
        sub_final = time.time()
        print "  Triaxiality time: ", sub_final - sub_initial
                       
        print sh.coefficients
        if config.reconstruct:
            sub_initial = time.time()
            sh.reconstruct(config.bins)
            sub_final = time.time()
            print "  Reconstruction time: ", sub_final - sub_initial
            
            sub_initial = time.time()
            sh.DensityPlotCompare(config.smoothing_scale, config.resolution)
            sub_final = time.time()
            print "  Density Plotting time: ", sub_final - sub_initial 
    
            sub_initial = time.time()
            sh.mserror()
            sub_final = time.time()
            print "  Mean Square Error time: ", sub_final - sub_initial
            
            sub_initial = time.time()
            sh.PeakSN()
            sub_final = time.time()
            print "  Peak Signal to Noise time: ", sub_final - sub_initial   
     
            print sh.rebuilt
            
            sub_initial = time.time()
            sh.Compare()
            sub_final = time.time()
            print "  Comparison time: ", sub_final - sub_initial   
                              
    final = time.time()
    print "Overall Time Taken: ", final - initial