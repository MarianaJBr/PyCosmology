'''
Created on May 7, 2012

@author: smurray
'''
import numpy as np

from PyCosmology.sims.fort.read_sim import simops as fsims
from PyCosmology.sims.fort.groupstructure import groupstructure as gs
import plotting

import os, errno
import matplotlib.pyplot as plt
class sims(object):
    """
    A class to import simulation data and group data from gadget runs. Basic plotting/manipulation.
    """
    
    def __init__(self,filename,group_filename=None,ids_filename=None,n_groups=1,reverse=False,centring_mode=None,loud=False,subgroup_filename=None):
        """
        Imports the simulation and detects groups and subgroups if specified.
        """
        #Set the FORTRAN module parameter whether to write out extra information.
        fsims.loud = loud
        
        # readsim reads in the simulation and defines fsims.pos and fsims.vel amongst others.
        fsims.readsim(filename)

        # find_groups locates groups in the simulation and saves the array fsmims.group_number,
        # which contains the number of the group for each particle.  
        if group_filename:
            fsims.find_groups(group_filename,ids_filename)
            
            # Modify the number of groups to work with in the stats part.
            if n_groups is 'all':
                self.n_groups = fsims.ngroups
            else:
                self.n_groups = n_groups
             
            #If a subgroup filename was given, read in all the subgroups, saving a flagged subgroup number to each particle.   
            if subgroup_filename:
                fsims.read_subgroupv(subgroup_filename)
                self.subgroups_pos = []
            
            #Initialize all variables for later use.    
            self.groups_pos = []
            self.groups_vel = []
            self.groups_pos_f = []
            self.groups_vel_f = []
            self.axis_b = []
            self.axis_c = []
            self.axis_d = []
            self.original_group_centres = []
            
            for i in range(self.n_groups):
                #Save the particles to their respective groups
                group_pos, group_vel = fsims.specify_groups(i+1,fsims.grouplen[i])
                self.groups_vel = self.groups_vel + [np.asfortranarray(group_vel)]
                self.groups_pos = self.groups_pos + [np.asfortranarray(group_pos)]
                
                #Save the original group centre (the first particle position of each group)
                self.original_group_centres = self.original_group_centres + [np.array(group_pos[:,0])]

                # Recentre each group around its first particle
                fsims.centre(self.groups_pos[i],centring_mode)
                
                #Rotate each group so its longest axis is the x-axis.
                fsims.rotate(self.groups_pos[i],np.asfortranarray(self.groups_vel[i]),[0.0,0.0,0.0],0.0)

                #Save axis ratios for further use.
                self.axis_b.append(float(fsims.axis_ratio_b))
                self.axis_c.append(float(fsims.axis_ratio_c))
                self.axis_d.append(float(fsims.axis_ratio_d))
                 
                #Locate subgroups.Particles are in no particular order though!
                self.subgroups_pos = self.subgroups_pos + [[]]
                
                #Arrange the subgroup numbers contained in each group in descending order so that non-grouped particles are last.
                indices = list(set(fsims.subgroup_number[fsims.group_number == i+1]))
                indices.sort(reverse=True)
                for index in indices:
                    self.subgroups_pos[i] = self.subgroups_pos[i] + [np.array(self.groups_pos[i])[:,fsims.subgroup_number[fsims.group_number == i+1]== index]]


        #Make the directory structure that files will be saved to.        
        self.make_directories(filename)
        
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
        
        self.group_density_dir = "GroupData/DensityPlots"
        try:
            os.makedirs(self.group_density_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
            
        self.group_profile_dir = "GroupData/DensityProfiles"
        try:
            os.makedirs(self.group_profile_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
        self.group_velocity_dir = "GroupData/VelocityDispersion"
        try:
            os.makedirs(self.group_velocity_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise        
            
        self.group_stats_dir = "GroupData/Statistics"
        try:
            os.makedirs(self.group_stats_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise  
            
        self.group_structure_shells_dir = "GroupData/StructureInShells"
        try:
            os.makedirs(self.group_structure_shells_dir)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise  
        self.nonparam = "GroupData/StructureInShellsNonParam"
        try:
            os.makedirs(self.nonparam)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise 
                       
    def plot_group_pos(self):
        """
        Convenience function for plotting the group positions as density maps all at once.
        """
        print "Plotting Group Positions"
        for i,group in enumerate(self.groups_pos):
            #If there are subgroups, we plot the subgroups as overlays with different colours.
            if self.subgroups_pos:
                plotting.SpatialPlot(group[0,:],group[1,:],self.group_density_dir+'/group_'+str(i),smoothing='dot',subgroups=self.subgroups_pos[i])
            else:
                plotting.SpatialPlot(group[0,:],group[1,:],self.group_density_dir+'/group_'+str(i),smoothing='dot')
                
    def plot_group_density_profile(self):
        """
        Convenience function for plotting the density profiles for all maintained groups.
        """
        print "Plotting Group Density Profiles"
        for i,group in enumerate(self.groups_pos):
            dens,edges = fsims.density_profile(group,10,0.01)
            plotting.DensityProfile(edges[:-1],dens,self.group_profile_dir+'/group_'+str(i))
        
    def plot_velocity_disp(self):
        """
        Convenience function for plotting the velocity dispersion profiles for all maintained groups.
        """
        print "Plotting Velocity Dispersions"
        for i,group in enumerate(self.groups_vel):
            plotting.VelocityDispersion(group.T, self.group_velocity_dir+'/group_'+str(i))


    def AxisRatioDistribution(self):
        """
        Makes a scatter plot of the axis ratios of all groups maintained.
        """
        print "Plotting Distribution of Axis Ratios"

        plt.clf()
        fig = plt.figure(1, figsize=(5.5,5.5))
        plt.suptitle("Axis ratio scatter for the groups (N = "+str(len(self.axis_b))+')')

        plt.xlabel("Second Longest to Longest Axis")
        plt.ylabel("Shortest to Second Longest Axis")
        from mpl_toolkits.axes_grid1 import make_axes_locatable
    
        # the scatter plot:
        axScatter = plt.subplot(111)
        axScatter.scatter(self.axis_b, self.axis_d)
        
        # create new axes on the right and on the top of the current axes
        # The first argument of the new_vertical(new_horizontal) method is
        # the height (width) of the axes to be created in inches.
        divider = make_axes_locatable(axScatter)
        axHistx = divider.append_axes("top", 1.2, pad=0.1, sharex=axScatter)
        axHisty = divider.append_axes("right", 1.2, pad=0.1, sharey=axScatter)
        
        # make some labels invisible
        plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
                 visible=False)
        
        # now determine nice limits by hand:
        #binwidth = 0.25
        #xymax = np.max( [np.max(np.fabs(self.axis_b)), np.max(np.fabs(self.axis_d))] )
        #lim = ( int(xymax/binwidth) + 1) * binwidth
        
        #bins = np.arange(-lim, lim + binwidth, binwidth)
        axHistx.hist(self.axis_b, bins=20)
        axHisty.hist(self.axis_d, bins=20, orientation='horizontal')
        
        # the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
        # thus there is no need to manually adjust the xlim and ylim of these
        # axis.
        
        #axHistx.axis["bottom"].major_ticklabels.set_visible(False)
        for tl in axHistx.get_xticklabels():
            tl.set_visible(False)
        axHistx.set_yticks([])
        
        #axHisty.axis["left"].major_ticklabels.set_visible(False)
        for tl in axHisty.get_yticklabels():
            tl.set_visible(False)
        axHisty.set_xticks([])

        plt.draw()
        plt.show()
        plt.savefig(self.group_stats_dir+'/AxisRatioScatter')
        
        plt.clf()
        plt.figure()
        plt.suptitle("Axis ratio distributions for "+str(len(self.axis_b))+" groups.")
        plt.subplot(311)
        plt.hist(self.axis_b,bins=20)
        plt.xlabel("b/a")
        plt.subplot(312)
        plt.hist(self.axis_c,bins=20)
        plt.xlabel("c/a")
        plt.subplot(313)
        plt.hist(self.axis_d,bins=20)
        plt.xlabel("d/c")
        plt.savefig(self.group_stats_dir+'/AxisRatioDistributions')
        
        
    def Centre_of_mass_offsets(self):
        """
        Calculates the distribution of centre of mass offsets (not right yet)
        """
        print "Plotting Centre of Mass Offsets"

        fof_offset = []
        sphere_offset = []   
        for i,group in enumerate(self.groups_pos):
            #Don't worry about the last group, which will be the ungrouped particles.
            if i == len(self.groups_pos)-1:
                continue
            
            fof_offset = fof_offset + [fsims.centre_of_mass_offset_fof(group)]
            sphere_offset = sphere_offset + [fsims.centre_of_mass_offset_sphere(np.asfortranarray(self.original_group_centres[i]))]
        plt.clf()
        plt.title("Centre of Mass Offsets For FoF and Spherical Overdensity Groups")
        plt.xlabel("FoF Mass Offset, log(Mpc)")
        plt.ylabel("Spherical Overdensity Offset, log(Mpc)")
        plt.plot(np.log(fof_offset),np.log(sphere_offset),'bo')
        plt.savefig(self.group_stats_dir+"/Centre_of_Mass_Offsets.png")
        
    def structure_in_shells(self, bins):
        """
        Calculates the ratio of subgroup particles to non-subgroup particles for a shell.
        """
        print "Calculating Ratio of Subgroup Particles to Background Particles in Subshells of Groups"
        
        
        average = np.zeros(bins-1)
        for k,group in enumerate(self.subgroups_pos):
            #Set up the radii
            end =  (3.0*fsims.grouplen[k]*fsims.massarr[1]/(4.0*np.pi*200.0*27.755))**(1./3.)
            radii = np.linspace(0.0, end, bins)
        
            ratio = []
            rad = []
            for i,radius in enumerate(radii):
                #We have n-1 bins where n is the bin edge number.
                if i==len(radii)-1:
                    continue
                
                n_subgroup_particles = 0
                for j,subgroup in enumerate(group):
                    #Subgroup particles are in the first n-1 subgroups
                    if j != len(group)-1:
                        n_subgroup_particles = n_subgroup_particles + fsims.number_in_shell(radius, radii[i+1],subgroup)
                    #Background particles are in the last 'subgroup'
                    else:
                        n_background_particles = fsims.number_in_shell(radius, radii[i+1],subgroup)
                        
                if n_subgroup_particles == 0  and n_background_particles == 0:
                    ratio = ratio + [0]
                else:
                    ratio = ratio + [float(n_subgroup_particles)/(float(n_background_particles)+float(n_subgroup_particles))]
                rad = rad + [(radius+radii[i+1])/(2.0*end)]
            average = average + fsims.grouplen[k]*np.array(ratio)
            
            plt.clf()
            plt.title("Substructure in Shells")
            plt.xlabel("Radius")
            plt.ylabel("Structure Fraction")
            plt.plot(rad,ratio)
            plt.savefig(self.group_structure_shells_dir+'/group_'+str(k))
        
        #Also plot an average, weighted by number of particles in groups.
        average = average/(self.n_groups*np.sum(fsims.grouplen))
        plt.clf()
        plt.title("Substructure in Shells On Average")
        plt.xlabel("Radius - fraction of Virial Radius")
        plt.ylabel("Structure Fraction")
        plt.plot(rad,average)
        plt.savefig(self.group_structure_shells_dir+'/all_groups')
        
    def non_parametric_shell_structure(self,n_radial_bins,n_angles,nside):
        """
        Does a non-parametric test for structure in radial shells of groups and clusters.
        """
        print "Performing Non-Parametric Variance-Based Test For Structure In Groups"
        
        for i,group in enumerate(self.groups_pos):
            variances, poisson, rads_c, rads, angles_c = gs.structure_test(group,n_radial_bins, n_angles,nside)

            plt.clf()
            plt.figure(figsize=(10.0,20.0))
            plt.subplots(nrows=3, ncols=1)
            plt.subplots_adjust(top=0.9, bottom = 0.1, hspace = 0.6)
            plt.suptitle("Varying Radius to Keep Constant Volume")
            
            plt.subplot(311)
            plt.title("Log Density Variance in Radial Shells on Varying Scales ")
            plt.imshow(np.log10(variances[1,:,:]),interpolation='Gaussian', aspect = float(n_angles)/float(n_radial_bins) )
            plt.xlabel("Angle Sizes")
            plt.ylabel("Radial Bins (Varied non-linearly)")
            plt.colorbar()

            plt.subplot(313)
            plt.title("Density Variance on Varying Scales Averaged Over Shells")
            shell_av_var = np.sum(variances[1,:,:],axis=0)/n_radial_bins
            shell_av_pos = np.sum(poisson[1,:,:],axis=0)/n_radial_bins
            plt.plot(angles_c,np.log10(shell_av_var), 'b-')
            plt.plot(angles_c,np.log10(shell_av_pos), 'g-')
            plt.xlabel("Angle Sizes")
            plt.ylabel("Log of Density Variance")
            for rbin in range(n_radial_bins):
                plt.plot(angles_c,np.log10(variances[1,rbin,:]),'m-')

            plt.subplot(312)
            plt.title("Density Variance in Radial Shells Average Over Scales")
            scale_av_var = np.sum(variances[1,:,:],axis=1)/n_angles
            scale_av_pos = np.sum(poisson[1,:,:],axis=1)/n_angles
            plt.plot(rads,np.log10(scale_av_var),'b-')
            plt.plot(rads,np.log10(scale_av_pos),'g-')
            plt.xlabel("Centre of Radial Shell as Fraction of Virial Radius")
            plt.ylabel("Log of Density Variance")
            for angle in range(n_angles):
                plt.plot(rads,np.log10(variances[1,:,angle]),'m-')
                
            plt.savefig(self.nonparam+"/group_"+str(i)+'_varRadius')
            
            plt.clf()
            plt.figure(figsize=(10.0,20.0))
            plt.subplots(nrows=2, ncols=1)
            plt.subplots_adjust(top=0.9, bottom = 0.1, hspace = 0.6)
            plt.suptitle("Varying Angles to Keep Constant Volume")
            
            plt.subplot(211)
            plt.title("Log Density Variance in Radial Shells on Varying Scales ")
            plt.imshow(np.log10(variances[0,:,:]),interpolation='Gaussian', aspect = float(n_angles)/float(n_radial_bins))
            plt.xlabel("Angle Sizes (Varied Non-Linearly")
            plt.ylabel("Radial Bins")
            plt.colorbar()
            
            plt.subplot(212)
            plt.title("Density Variance in Radial Shells Average Over Scales")
            scale_av_var = np.sum(variances[0,:,:],axis=1)/n_angles
            scale_av_pos = np.sum(poisson[0,:,:],axis=1)/n_angles
            plt.plot(rads_c,np.log10(scale_av_var),'b-')
            plt.plot(rads_c,np.log10(scale_av_pos),'g-')
            plt.xlabel("Centre of Radial Shell as Fraction of Virial Radius")
            plt.ylabel("Log of Density Variance")
            for angle in range(n_angles):
                plt.plot(rads_c,np.log10(variances[0,:,angle]),'m-')
                
            plt.savefig(self.nonparam+"/group_"+str(i)+'_varAngle')
            