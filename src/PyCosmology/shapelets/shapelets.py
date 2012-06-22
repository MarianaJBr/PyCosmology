'''
Created on Jun 5, 2012

@author: Steven
'''

from PyCosmology.sims.fort.read_sim import simops
import PyCosmology.sims.plotting as plotting
from PyCosmology.shapelets.fort.shapelets import shapelets as fc
import numpy as np
from scipy.misc import comb
import matplotlib.pyplot as plt
import os, errno
from bigfloat import erf
#import matplotlib
#matplotlib.use("Qt4Agg")
#import matplotlib.pyplot as plt


class ShapeletOperations_DMH(object):
    """
    This class contains methods for evaluating shapelet co-efficients and plotting
    results based on this.
    """
    
    def __init__(self,sim_file,groups_file,group_ids_file, n=1,reverse=False,test=False,basis_test=False,n_test=51):
        """
        Read in a sim-file for use with the shapelet analysis.
        
        PARAMETERS:
        sim_file: file containing the simulation snapshot
        groups_file: file containing the location of groups
        group_ids_file: file containing the ids of the group members
        
        n: optional, default 1, number of groups to retain for further analysis (from largest to smallest)
        reverse: optional, default False, whether to use the n SMALLEST groups rather than largest.
        """
        
        if basis_test:
            self.test_bases(n_test)
            return
            
        self.test = test
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
            
        self.make_directories(sim_file)

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
   
        os.chdir(self.sim_dir)
        
        self.group_dir = "GroupData"
        try:
            os.makedirs(self.group_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
        
        self.shapelet_dir = "GroupData/Shapelets"
        try:
            os.makedirs(self.shapelet_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
        os.chdir(self.shapelet_dir)
        
        self.reconstruct_compare_dir = "ReconstructionComparison"
        try:
            os.makedirs(self.reconstruct_compare_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
          
 
    def find_coeffs(self,bins):
        """
        Finds the coefficients for the shapelets
        """
        
        print "FINDING SHAPELET COEFFICIENTS"
        self.coefficients = []
        self.beta = []
        self.nmax = (bins-3)/2
        self.hist = []
        self.x_max = []
        if not self.test:
            for group in self.groups_pos:
                x_max = 2.0*np.max(group[0,:])
                self.x_max = self.x_max + [x_max]
    
                hist,edges = np.histogramdd(group.T,bins=bins,range=[[-x_max,x_max],[-x_max,x_max],[-x_max,x_max]])
                del edges
    
                hist = np.log10(1.0+simops.massarr[1]*hist)
                self.hist= self.hist + [hist]
    
                self.coefficients = self.coefficients + [np.asfortranarray(fc.coeff_cube(np.asfortranarray(hist),x_max.astype(np.float128)))]
                self.beta = self.beta+[x_max/np.sqrt(2*bins)]
                
        else:
            self.test_sphere()
            x_max = 2.0
            self.x_max = self.x_max + [x_max]
    
            hist,edges = np.histogramdd(self.test_group,bins=bins,range=[[-x_max,x_max],[-x_max,x_max],[-x_max,x_max]])
            del edges
    
            hist = np.log10(1.0+simops.massarr[1]*hist)
            self.hist= self.hist + [hist]
    
            self.coefficients = self.coefficients + [np.asfortranarray(fc.coeff_cube(np.asfortranarray(hist),x_max))]
            self.beta = self.beta+[fc.beta]            

    def reconstruct(self,bins):
        """
        Reconstructs the 'original' data using the basis functions and amplitudes
        """
        
        print "RECONSCTRUCTING GROUP DENSITY"
        self.rebuilt = []
        for i,coefficients in enumerate(self.coefficients):
            self.rebuilt = self.rebuilt + [fc.reconstruct_cube(bins,coefficients,(-self.x_max[i],self.x_max[i]))]
            
                            
    def zeroth_moment(self):
        """
        Finds the zeroth moment of the shape based on the coefficients found. Associated with total mass.
        """
        
        print "CALCULATING ZEROTH MOMENT"
        self.M_0 = []
        a = np.arange(0,self.nmax,2)

        
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
        
        """
        
        print "CALCULATING CENTROID"
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
        
            self.centroid = self.centroid + [(x1,y1,z1)]
            print self.centroid[i]
        
    def rmsRadius(self):
        """
        Calculates the rms radius defined in Fluke. et. al.
        """
        print "Calculating RMS Radius Value"
        a = np.arange(0,self.nmax,2)
        
        self.rms = []
        for i,beta in enumerate(self.beta):
            rms = 0.0
            for n1 in a:
                for n2 in a:
                    for n3 in a:
                        rms = rms + self.coefficients[i][n1,n2,n3]*(n1+n2+n3+1.5)*self.U(n1+n2+n3)*self.W(n1,n2,n3)
                        
            rms = (2.0*np.pi**(3./4.)*beta**(3.5)/self.M_0[i]) * rms
            self.rms = self.rms + [rms]  
            
            print rms    
            
    def InertiaTensor(self):
        """
        Calculate the moment of inertia via Fluke. et. al.
        """
        
        print "Calculating the Moment of Inertia Tensor via shapelet coefficients."
        self.InertiaTens = []
        a = np.arange(0,self.nmax,2)
        b = np.arange(1,self.nmax,2)
        
        for i,beta in enumerate(self.beta):
            tensor = np.zeros((3,3))
            
            for n1 in a:
                for n2 in a:
                    for n3 in a:
                        tensor[0,0] =  tensor[0,0] + self.coefficients[i][n1,n2,n3]*(n2+n3+1)*self.U(n1+n2+n3)*self.W(n1,n2,n3)
                        tensor[1,1] =  tensor[1,1] + self.coefficients[i][n1,n2,n3]*(n1+n3+1)*self.U(n1+n2+n3)*self.W(n1,n2,n3)
                        tensor[2,2] =  tensor[2,2] + self.coefficients[i][n1,n2,n3]*(n2+n1+1)*self.U(n1+n2+n3)*self.W(n1,n2,n3)
            
 
                
            for n3 in a:
                for n1 in b:
                    for n2 in b:
                        tensor[0,1] = tensor[0,1] -  self.coefficients[i][n1,n2,n3]*np.sqrt((n1+1)*(n2+1))*self.U(n1+n2+n3)*self.W(n1+1,n2+1,n3)
                        tensor[0,2] = tensor[0,2] -  self.coefficients[i][n1,n3,n2]*np.sqrt((n1+1)*(n2+1))*self.U(n1+n2+n3)*self.W(n1+1,n2+1,n3)
                        tensor[1,2] = tensor[1,2] -  self.coefficients[i][n3,n2,n1]*np.sqrt((n1+1)*(n2+1))*self.U(n1+n2+n3)*self.W(n1+1,n2+1,n3)
                      
            tensor[1,0] = tensor[0,1]
            tensor[2,0] = tensor[0,2]
            tensor[2,1] = tensor[1,2]
                      
            tensor = tensor* 2.0*np.pi**(3./4.) * beta**3.5
            
            self.InertiaTens = self.InertiaTens + [tensor]  
            
            print tensor
    def Triaxiality(self):
        """
        Calculate the triaxiality from the appendix of Fluke
        """
        self.triax_shape = []
        for inertia in self.InertiaTens:
            self.triax_shape = self.triax_shape +[(inertia[1,1]-inertia[0,0])/(inertia[2,2]-inertia[0,0])]
            
        print self.triax_shape
    def PeakSN(self):
        """
        Calculate the peak signal to noise ratio for each group
        """
        self.Ps = []
            
        for i,cube in enumerate(self.hist):
            self.Ps = self.Ps + [20*np.log10(np.max(cube)/np.sqrt(self.MSError[i]))]
            
            print "Peak Signal to noise: ", self.Ps[i]  
            
    def mserror(self):
        """
        calculate the mean square error in the reconstrution.
        """
        self.MSError = []
        
        for i,cube in enumerate(self.hist):
            ms = 0.0
            for j in range(cube.shape[0]):
                for k in range(cube.shape[0]):
                    for l in range(cube.shape[0]):
                        ms = ms + (cube[j,k,l] - self.rebuilt[i][j,k,l])**2
            ms = ms/cube.shape[0]**3
            
            self.MSError = self.MSError + [ms]
            
            print "Mean square error: ", ms
    def Compare(self):
        """
        Calculates quantities that compare the original to reconstructed groups
        """
        print "Doing overall comparison of reconstruction."
        self.compare = []
        for i,cube in enumerate(self.hist):
            sum_I = np.sum(cube)
            sum_S = np.sum(self.rebuilt[i])
            compare = sum_S/sum_I -1.0
            
            self.compare = self.compare+[compare]
            print "Summed quantities difference: ", compare    
                                  
    def DensityPlotCompare(self,smoothing_scale, resolution):
        """
        Plot the density projection of the original groups and their reconstructions to compare.
        """
        print "PLOTTING DENSITY MAPS FOR COMPARISON"
        if not self.test:
            for i,group in enumerate(self.groups_pos):
                plt.clf()
                plt.figure()
                plt.suptitle("Density Field for Group "+str(i))
                plt.subplot(121)
                plt.title("Original")
                #window_size = np.ceil(smoothing_scale*resolution)
                x = group[0,:]
                y = group[1,:]
                binned_data = np.log10(np.histogram2d(x,y, resolution,normed=True)[0]+1.0)
                re_ordered_data = np.transpose(np.fliplr(binned_data))
                #print re_ordered_data.shape
                #smoothed_data = plotting.sgolay2d(re_ordered_data,window_size=window_size,order=4)
                plt.imshow(re_ordered_data,extent=(min(x),max(x),min(y),max(y)),interpolation="Gaussian",aspect='equal')
                
                plt.subplot(122)
                plt.title("Reconstructed")
    
                density = np.sum(self.rebuilt[i],axis=2)
                #smoothed_data = plotting.sgolay2d(re_ordered_data,window_size=window_size,order=4)
                plt.imshow(density,interpolation="Gaussian",aspect='equal',extent=(-self.x_max[i],self.x_max[i],-self.x_max[i],self.x_max[i]))
                plt.show()
            
  
                plt.savefig(self.reconstruct_compare_dir+'/Group_'+str(i))
        else:
            plt.clf()
            plt.figure()
            plt.suptitle("Density Field for Test Group")
            plt.subplot(121)
            plt.title("Original")
                #window_size = np.ceil(smoothing_scale*resolution)
            x = np.array(self.test_group[0])
            y = np.array(self.test_group[1])
            binned_data = np.log10(np.histogram2d(x,y, resolution,normed=True)[0]+1.0)
            re_ordered_data = np.transpose(np.fliplr(binned_data))
                #print re_ordered_data.shape
                #smoothed_data = plotting.sgolay2d(re_ordered_data,window_size=window_size,order=4)
            plt.imshow(re_ordered_data,extent=(min(x),max(x),min(y),max(y)),interpolation="Gaussian",aspect='equal')
                
            plt.subplot(122)
            plt.title("Reconstructed")
    
            density = np.sum(self.rebuilt[0],axis=2)
                #smoothed_data = plotting.sgolay2d(re_ordered_data,window_size=window_size,order=4)
            plt.imshow(density,interpolation="Gaussian",aspect='equal',extent=(-self.x_max[0],self.x_max[0],-self.x_max[0],self.x_max[0]))                    
            plt.savefig(self.reconstruct_compare_dir+'/TEST_GROUP')

        
    def test_sphere(self):
        """
        Creates a mock group in the shape of a power-law sphere.
        """
        print "Making a test sphere"  
        self.test_group = np.zeros((3,5000))
        self.test_group[0,:] = np.random.rand(5000)
        self.test_group[1,:] = np.random.rand(5000)*2*np.pi
        self.test_group[2,:] = np.random.rand(5000)*np.pi
            
        self.test_group = self.Convert_to_Cartesian(self.test_group)
        
    def Convert_to_Cartesian(self,pos):
        """
        Convert a set of spherical polar co-ordinates into cartesian 3-space
        """
  
        x_pos = pos[0,:]*np.cos(pos[2,:])*np.sin(pos[1,:])
        y_pos = pos[0,:]*np.sin(pos[2,:])*np.sin(pos[1,:])
        z_pos = pos[0,:]*np.cos(pos[2,:])
            
        return [x_pos,y_pos,z_pos]
    
    def test_bases(self,n_test):
        """
        Use the fortran code to test basis vector outputs
        """
        
        cube0,cube1,cube2,cube3 = fc.test(n_test)
        
        plt.clf()
        plt.title("Basis Vector (0,0,0)")
        plt.imshow(cube0[:,:,n_test/2], interpolation='bilinear', origin='lower',
                        extent=(-3,3,-3,3))
        plt.contour(cube0[:,:,n_test/2],origin='lower',linewidths=2,extent=(-3,3,-3,3))
        plt.savefig('/Users/Steven/Documents/cube_0_test') 
        
        plt.clf()
        plt.title("Basis Vector (0,1,0)")
        plt.imshow(cube1[:,:,n_test/2], interpolation='bilinear', origin='lower',
                        extent=(-3,3,-3,3))
        plt.contour(cube1[:,:,n_test/2],origin='lower',linewidths=2,extent=(-3,3,-3,3))
        plt.savefig('/Users/Steven/Documents/cube_1_test') 
        
        plt.clf()
        plt.title("Basis Vector (2,0,2)")
        plt.imshow(cube2[:,:,n_test/2], interpolation='bilinear', origin='lower',
                        extent=(-3,3,-3,3))
        plt.contour(cube2[:,:,n_test/2],origin='lower',linewidths=2,extent=(-3,3,-3,3))
        plt.savefig('/Users/Steven/Documents/cube_2_test')  
        
        plt.clf()
        plt.title("Basis Vector (1,2,4)")
        plt.imshow(cube3[:,:,n_test/2], interpolation='bilinear', origin='lower',
                        extent=(-3,3,-3,3))
        plt.contour(cube3[:,:,n_test/2],origin='lower',linewidths=2,extent=(-3,3,-3,3))
        plt.savefig('/Users/Steven/Documents/cube_3_test')   