'''
Created on May 8, 2012

@author: smurray
'''
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import numpy as np

import Smoothing

def SpatialPlot(x,y,filename,smoothing="dot",condition ="",smoothing_scale=0.05,resolution=300):
    """
    Plots the spatial positions of objects in two-dimensions, using smoothing if dictated.
    """
    plt.clf()
    if smoothing is"dot":
        plt.title("Object Positions "+condition)
        plt.plot(x,y,'.',linewidth=0.01,markersize=0.01)
        plt.axes().set_aspect('equal', 'datalim')
        
    elif smoothing is "sgolay":
        plt.title("Density Field "+condition)
        window_size = np.ceil(smoothing_scale*resolution)
        
        binned_data = np.log(np.histogram2d(x, y, resolution,normed=True)[0]+1E-8)
        re_ordered_data = np.transpose(np.fliplr(binned_data))
        smoothed_data = Smoothing.sgolay2d(re_ordered_data,window_size=window_size,order=4)
        plt.imshow(smoothed_data,extent=(min(x),max(x),min(y),max(y)))
        plt.axes().set_aspect('equal', 'datalim')
        plt.show()
        
    else:
        plt.title("Density Field "+condition)
        
        binned_data = np.log(np.histogram2d(x, y, resolution,normed=True)[0]+1E-8)
        re_ordered_data = np.transpose(np.fliplr(binned_data))
        plt.imshow(re_ordered_data,extent=(min(x),max(x),min(y),max(y)),interpolation=smoothing)
        plt.axes().set_aspect('equal', 'datalim')
        plt.show()
        

    
    plt.savefig(filename)
    
    
def DensityProfile(r,dens,filename):
    """
    Outputs a plot of the density profile of a sample
    """
    
    plt.clf()
    plt.title("Density profile of the halo")
    plt.plot(r,np.log(dens),'ro',label="Measured Density")
    #plt.plot(r_NFW,dens_NFW,'b-', label="NFW Best Fit")
    plt.legend()
    plt.savefig(filename)
    
def VelocityDispersion(vel,filename):
    """
    Plots a histogram of velocities on given axes
    """
    
    plt.clf()
    plt.title("Velocity Dispersion")
    plt.hist(vel, 50,histtype='step',label=["Axis 1", "Axis 2", "Axis 3"], color=["r","b",'g'] )
    plt.legend()
    plt.savefig(filename)