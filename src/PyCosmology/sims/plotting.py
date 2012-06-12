'''
Created on May 8, 2012

@author: smurray
'''
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import numpy as np

import scipy.signal



def sgolay2d ( z, window_size, order, derivative=None):
    """
    """
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2)  / 2.0

    if  window_size % 2 == 0:
        print 'window_size must be odd, adding one.'
        window_size= window_size+1

    if window_size**2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial. 
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ... 
    # this line gives a list of two item tuple. Each tuple contains 
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]

    # coordinates of points
    ind = np.arange(-half_size, half_size+1, dtype=np.float64)
    dx = np.repeat( ind, window_size )
    dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )

    # build matrix of system of equation
    A = np.empty( (window_size**2, len(exps)) )
    for i, exp in enumerate( exps ):
        A[:,i] = (dx**exp[0]) * (dy**exp[1])

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    Z = np.zeros( (new_shape) )
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] =  band -  np.abs( np.flipud( z[1:half_size+1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band  + np.abs( np.flipud( z[-half_size-1:-1, :] )  -band )
    # left band
    band = np.tile( z[:,0].reshape(-1,1), [1,half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs( np.fliplr( z[:, 1:half_size+1] ) - band )
    # right band
    band = np.tile( z[:,-1].reshape(-1,1), [1,half_size] )
    Z[half_size:-half_size, -half_size:] =  band + np.abs( np.fliplr( z[:, -half_size-1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z

    # top left corner
    band = z[0,0]
    Z[:half_size,:half_size] = band - np.abs( np.flipud(np.fliplr(z[1:half_size+1,1:half_size+1]) ) - band )
    # bottom right corner
    band = z[-1,-1]
    Z[-half_size:,-half_size:] = band + np.abs( np.flipud(np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band )

    # top right corner
    band = Z[half_size,-half_size:]
    Z[:half_size,-half_size:] = band - np.abs( np.flipud(Z[half_size+1:2*half_size+1,-half_size:]) - band )
    # bottom left corner
    band = Z[-half_size:,half_size].reshape(-1,1)
    Z[-half_size:,:half_size] = band - np.abs( np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band )

    # solve system and convolve
    if derivative == None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -c, mode='valid')
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid')
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid'), scipy.signal.fftconvolve(Z, -c, mode='valid')

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
        
        binned_data = np.log10(np.histogram2d(x, y, resolution,normed=True)[0]+1.0)
        re_ordered_data = np.transpose(np.fliplr(binned_data))
        smoothed_data = sgolay2d(re_ordered_data,window_size=window_size,order=4)
        plt.imshow(smoothed_data,extent=(min(x),max(x),min(y),max(y)))
        plt.axes().set_aspect('equal', 'datalim')
        plt.show()
        
    else:
        plt.title("Density Field "+condition)
        
        binned_data = np.log10(np.histogram2d(x, y, resolution,normed=True)[0]+1.0)
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
    plt.plot(r,np.log10(1.0+dens),'ro',label="Measured Density")
    #plt.plot(r_NFW,dens_NFW,'b-', label="NFW Best Fit")
    plt.xlabel("Radius from the centre of the group")
    plt.ylabel("log_10(1+rho)")
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