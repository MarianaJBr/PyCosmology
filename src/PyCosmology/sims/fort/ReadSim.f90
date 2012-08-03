include 'Auxiliary.f90'

module SimOps
	use AuxiliaryOperations
	use pix_tools !Have include directory and library directory in compile options.
	implicit none

	!Fundamental Parameters. 
	real(8), parameter		:: pi = 3.141592653589	!Fundamental Constant.
	real(8), parameter		:: delta_vir = 200.d0	!Critical Overdensity defining virial radius.
	real(8), parameter		:: rho_c = 27.75d0		!Critical Density of the Universe.
	
	!Important Properties of the Simulation
	real, allocatable		:: pos(:,:)			!3D Positions of all particles
	real, allocatable		:: vel(:,:)			!3D Velocities of all particles
	integer  				:: FlagSfr			!
	integer					:: FlagFeedback		!
	integer					:: Nall(6)			!
	integer					:: FlagAge			! 
	integer					:: FlagCooling 		!
	integer  				:: FlagMetals		!
	integer					:: NallHW(6)		!
	integer					:: FlagEntrICs		!
	integer					:: npart(6)			!Number of various types of particles. DM only will have npart(2) /=0
	integer					:: NumFiles			!Number of simulation files.
	real(8)  				:: massarr(6)		!Mass of each particle type. DM only has massarr(2) only
	real(8)					:: time				!
	real(8)					:: BoxSize			!Box Size of the Simulation
	real(8)					:: OmegaLambda		!Cosmological Constant of Simulation.
	real(8)	 				:: HubbleParam		!Hubble Parameter of Simultion
	real(8)					:: Omega0			!Fraction of Criticial Density of Current Universe
	real(8)					:: Redshift			!Redshift of the Simulation.
	
	!Important Group Properties
	integer					:: ngroups			!Total number of groups in sim.
	integer, allocatable	:: group_number(:)	!The group number of each particle (in order of largest to smallest group)
	integer,allocatable		:: groupoffset(:)	!The offset of each group.
	integer,allocatable		:: grouplen(:)		!The size of each group
	integer, allocatable	:: ids(:)			!ID's of all particles IN GROUPS
	real					:: axis_ratio_b		!2ndAxis/1stAxis of a group
	real					:: axis_ratio_c		!3rdAxis/1stAxis of a group
	real					:: axis_ratio_d		!3rdAxis/2ndAxis of a group
	
	
	!Important Subgroup Properties.
    integer					:: nsub				!Number of subhalos.
    integer, allocatable	:: subparent(:)    !fofcat ID of parent halo. Allocated as nsub.
    integer, allocatable	:: subgroup_number(:)	!The subgroup number (ID) of every particle.
	
    !Run Parameters
	logical		:: loud = .FALSE.

	contains


  !----------------------------------------------------------------------
 	  ! ReadSim - Simply read in information from the simulation
  !----------------------------------------------------------------------		
    subroutine ReadSim(filename)    
    	character(len=*), intent(in)	:: filename	!Filename of the Simulation
    	
    	integer 						:: i		!Iterator.
    	integer, allocatable			:: ID(:)	!ID's of all particles.
    	
		open(unit=10,file=filename,form='UNFORMATTED')
        
		! Read Overall Simulation Properties. All saved in case needed.
		read(10) npart, massarr, time, Redshift, FlagSfr, FlagFeedback, Nall,&
				  & FlagCooling, NumFiles, BoxSize, Omega0, OmegaLambda,&
				  &HubbleParam,FlagAge, FlagMetals, NallHW, FlagEntrICs
		
		! Write interesting information out.
		write(*,*) "Number of particles = ", npart(2)
		write(*,*) "Time of snapshot: ", time
		write(*,*) "Redshift of snapshot: ", Redshift
		write(*,*) "BoxSize of Sim: ", BoxSize
		
		!Make sure the relevant arrays are not previously allocated
		! which may happen if reading in multiple sim files.
		if (allocated(pos))then
			deallocate(pos)
		end if
		if (allocated(vel))then
			deallocate(vel)
		end if
		
		!Allocate the appropriate arrays.
		allocate(pos(3,npart(2)),vel(3,npart(2)),ID(npart(2)))
		
		!Skip ahead to the ID's for automatic ordering.
		read(10)
		read(10)
				
		! Read in Simulation File values (only positions, velocities and ID's)
		read(10) (ID(i), i=1,npart(2))
		write(*,*) "Read Particle ID's"
		
		!Rewind to the positions and velocities.
		rewind(10)
		read(10)
		read(10) (pos(1,ID(i)),pos(2,ID(i)),pos(3,ID(i)), i=1,npart(2))
		write(*,*) "Read Positions"
		read(10) (vel(1,ID(i)),vel(2,ID(i)),vel(3,ID(i)), i=1,npart(2))
		write(*,*) "Read Velocities"
		close(10)
      
		!The ID's are not needed for anything else, so they can be deallocated.
		deallocate(ID)
    end subroutine
    
  !----------------------------------------------------------------------
 	  ! find_groups - establish which objects belong to which groups.
  !----------------------------------------------------------------------
    subroutine find_groups(filename,ids_file)  
		character(len=*), intent(in)	:: filename		!Group catalogue file
		character(len=*), intent(in)	:: ids_file 	!Particle ID's of all particles in groups.


    	integer 			:: k,i				!Iterators
    	integer				:: ntot				!Total number of particles in groups
    	integer				:: group_offset		!The offset from 0 of a particular group in the ID's file.
    	integer				:: offset_pos(1)	!The identifier of the next largest group
    	integer				:: group_ID			!The ID of a particular member of a group.

    	
    	!Open the FoF group catalogue
    	open(1,file=filename,status='unknown',form='unformatted')
    	          
    	!First get the total number of groups.
    	read(1) ngroups
    	
    	write(*,'(a,i6,a)') 'Reading ',ngroups,' groups...'
    	
    	!Allocate the main arrays.
    	if(allocated(grouplen))then
    		deallocate(grouplen)
    	end if
    	
    	allocate(grouplen(ngroups))
    	allocate(groupoffset(ngroups))
    	
    	!Read the information into the main arrays.
    	read(1) (grouplen(i),i=1,ngroups) 
    	write(*,*) 'Read group lengths'
    	read(1) (groupoffset(i),i=1,ngroups)
    	write(*,*) 'Read group offsets'
    	!read(1) (nsubgroups(i),i=1,ngroups) ! Number of subgroups -- ignore
    	!write(*,*) 'Read number of subgroups'
    	close(1)
    	
    	!Now Open the ID's file
    	open(1,file=ids_file,status='unknown',form='unformatted')
    	read(1) ntot          
    	write(*,*) "There are ", ntot, " particles in groups"
    	allocate(ids(ntot))
    	read(1) (ids(i),i=1,ntot)
    	close(1)    
    	
    	!=======================================================
    	!We now have  how big each group is, and where each group starts in the
    	!ID file, and the ID's of all the particles in the groups. We only need
    	!to flag each particle with its correct group now (actually building
    	!lists of groups will be left to the more flexible python).
    	!=======================================================
    		
    	!Allocate an array to store the group number (in descending size order)
    	!of each particle.
        allocate(group_number(npart(2)))
        
        !Initialize to zero, so particles with the zero flag are not in groups.
        group_number(:) = 0
    
        !The process here is to locate the largest group, and flag all particles
        !belonging to it as belonging to group 1, then find the next largest 
        !group and set them to 2 etc.
        !do k = 1,ngroups
        !	
        !	!Locate largest group
        !	offset_pos = maxloc(grouplen)
        !	group_offset = groupoffset(offset_pos(1))
        !  	
        !  	!Specify the particles of this group, returning an array of numbers
        !  	!corresponding to which group each is in.
        !  	do i=1,maxval(grouplen)
        !  	  group_ID = ids(group_offset+i)
        !  	  group_number(group_ID) = k
        !  	
        !  	end do
        !  	
        !  	!Set the current group size to negative so the next biggest group
        !  	!can be found simply.
        !  	grouplen(offset_pos(1)) = -1
        !end do
        	
        do k=1,ngroups
        	do i=1,grouplen(k)
        		group_ID = ids(groupoffset(k)+i)
        		group_number(group_ID) = k
        	end do
        end do
      end subroutine

  !----------------------------------------------------------------------
 	  ! specify_groups - determine the actual particle pos and vel for a group.
  !----------------------------------------------------------------------   
  	subroutine specify_groups(g_num,n,g_pos,g_vel)
  		integer, intent(in)		:: n 						!Size of group
  		integer, intent(in)		:: g_num					!The index of the group being specified
  		real, intent(out)		:: g_pos(3,n)				!The positions of the group
  		real, intent(out)		:: g_vel(3,n)				!The velocities of the group

  		integer					:: i,k						!Iterator
  		integer					:: group_ID					!The ID of a particular member of a group.
  		
  		do i = 1,n
        	group_ID = ids(groupoffset(g_num)+i)
        	g_pos(:,i) = pos(:,group_ID)
        	g_vel(:,i) = vel(:,group_ID)
  			!if(i==1)then
  			!	write(*,*) "First particle in group pos:", g_pos(:,i)
  			!end if
  		end do
  	end subroutine
  !----------------------------------------------------------------------
 	  ! specify_group_pos - determine the actual particle pos for a group.
  !----------------------------------------------------------------------   
  	subroutine specify_group_pos(g_num,n,g_pos)
  		integer, intent(in)		:: n 						!Size of group
  		integer, intent(in)		:: g_num					!The index of the group being specified
  		real, intent(out)		:: g_pos(3,n)				!The positions of the group

  		integer					:: i,k						!Iterator
  		
  		do i = 1,n
  			do k=1,npart(2)
  				if(group_number(k) == g_num)then
  					g_pos(:,i) = pos(:,k)
  					
  					exit
  				end if
  			end do
  			if(i==1)then
  				write(*,*) "First particle in group pos:", g_pos(:,i)
  			end if
  		end do
  	end subroutine
  !----------------------------------------------------------------------
 	  ! specify_group_vel - determine the actual particle vel for a group.
  !----------------------------------------------------------------------   
  	subroutine specify_group_vel(g_num,n,g_vel)
  		integer, intent(in)		:: n 						!Size of group
  		integer, intent(in)		:: g_num					!The index of the group being specified
  		real, intent(out)		:: g_vel(3,n)				!The velocities of the group
  		
  		integer					:: i,k						!Iterator
  		
  		do i = 1,n
  			do k=1,npart(2)
  				if(group_number(k) == g_num)then
  					g_vel(:,i) = vel(:,k)
  					
  					exit
  				end if
  			end do

  		end do
  	end subroutine	
  !----------------------------------------------------------------------
 	  ! centre - change co-ords of a group to its centre
  !----------------------------------------------------------------------        
        
  	subroutine centre(g_pos,n,mode)
  	!
  	!mode refers to where to place the centre (0,0,0).
  	!
  	!mode = 1 uses the first particle (centre of density)
  	!mode = 2 uses the average particle position (centre of mass)
  	!mode = 3 recentres with equal limits in both directions.
  	!	

  	
    	integer,intent(in)	:: n 					!The size of the group
    	integer,intent(in)	:: mode					!The mode to recentre in
  		real,intent(inout)	:: g_pos(3,n)			!The positions of the group.
  		
  		real				:: centre_of_density(3)	!The most bound particle (1st in the group)
  		real				:: c_o_m(3)				!Centre of Mass (average particle position)
  		real				:: equal_centre(3)		!Geometric Centre.
  		
  		integer				:: i					!Iterator
  		integer				:: temp_mode			!Temp mode fill in.

 	
  		
  		!Always centre about the first particle firstly, so that particles
  		!that are wrapped around the box don't affect the other measurements.
        centre_of_density = g_pos(:,1)
        
        ! Re-centre the group around the centre of density
        call ReCentre(g_pos, BoxSize,centre_of_density)
        if (loud) then
        write(*,"(A7,3es20.8,A1)") "CoD = (",&
                                    &centre_of_density,")"
        end if
        
        if (temp_mode == 2) then
        	
			c_o_m(:) = 0.000
			do i = 1,n
			  c_o_m(1) = c_o_m(1) + g_pos(1,i)
			  c_o_m(2) = c_o_m(2) + g_pos(2,i)
			  c_o_m(3) = c_o_m(3) + g_pos(3,i)
			end do
			c_o_m(:) = c_o_m(:)/real(n)
			if (loud)then
			write(*,"(A27,3es20.8,A1)") "Centre of Mass of group = (", c_o_m,")"
			end if
			! Re-centre the group around the centre
				call ReCentre(g_pos, BoxSize,c_o_m)
		end if
        
		if (temp_mode == 3) then
			equal_centre(1) = (maxval(g_pos(1,:))+minval(g_pos(1,:)))/2.d0
			equal_centre(2) = (maxval(g_pos(2,:))+minval(g_pos(2,:)))/2.d0
			equal_centre(3) = (maxval(g_pos(3,:))+minval(g_pos(3,:)))/2.d0
			! Re-centre the group around the centre
				call ReCentre(g_pos, BoxSize,equal_centre)
		end if
		
		if (loud) then
        write(*,*) "Group Re-centred"
        end if
     end subroutine
     
            
  
  !----------------------------------------------------------------------
 	  ! rotate - rotate data into moment of inertia defined coords plus angle
  !----------------------------------------------------------------------            
	subroutine rotate(g_pos, g_vel, n, u, theta)
		integer, intent(in)	::n				!The number of objects
		real, intent(inout)	::g_pos(3,n)	!The positions of the objects
		real, intent(inout)	::g_vel(3,n)	!The velocities of the objects
		real, intent(in)	::u(3)			!The vector around which to rotate (after rotation to the principle axes)
		real, intent(in)	::theta			!The angle by which to rotate (after rotation to the principle axes)
    
		
		integer	   			:: nrot			!A variable needed for the jacobi subroutine.
		integer				:: i			!Iterator
    	real      			:: IT(3,3)		!Moment of Inertia Tensor
    	real				:: eigval(3)	!Eigenvalues of the Inertia Tensor
		real      			:: eigvec(3,3)	!Eigenvectors of the Inertia Tensor
		real,allocatable	:: pos_r(:,:)	!Intermediate array for rotated co-ords.
		
        ! Calculate the Moment of Inertia tensor
        IT = Inertia(real(massarr(2)),g_pos)
        if (loud) then
        	write(*,*) "Moment of Inertia calculated"
        
        	write(*,*) "Moment of Inertia tensor:"
        	do i=1,3
            	write(*,*) IT(i,:)
            end do
        end if
         
        
        ! Diagonalize Inertia Tensor
          call jacobi(IT,3,3,eigval,eigvec,nrot)
          if (loud) then
          	  write(*,*) "Moment of Inertia tensor diagonalized"
          end if
          
        ! Re-order the eigs
          call EigOrder(eigval,eigvec)
          
        ! Return the ratios of axes
          axis_ratio_b = sqrt(eigval(2)/eigval(1))
          axis_ratio_c = sqrt(eigval(3)/eigval(1))
          axis_ratio_d = sqrt(eigval(3)/eigval(2))
          if (loud) then
          	  write(*,*) "Ratio of axes lengths b/a = ", axis_ratio_b
          	  write(*,*) "Ratio of axes lengths c/a = ", axis_ratio_c
          	  write(*,*) "Ratio of axes lengths c/b = ", axis_ratio_d
          end if
          
        ! Project data into co-ordinates defined by
            ! xbar = eigvec(1)*(x,y,z)
            ! ybar = eigvec(2)*(x,y,z)
            ! zbar = eigvec(3)*(x,y,z)
              
          do i=1,n
            g_pos(:,i) = matmul(eigvec,g_pos(:,i))
            g_vel(:,i) = matmul(eigvec,g_vel(:,i))
          end do
          
          allocate(pos_r(3,n))
          
          !Now the group has been rotated to the principle axes, and we can
          !perform another general rotation upon it.
          call GenRotate(g_pos,u,theta,pos_r)
          g_pos = pos_r
        
          call GenRotate(g_vel,u,theta,pos_r)
          g_vel = pos_r
          
      end subroutine
      
      
  !----------------------------------------------------------------------
 	  ! density_profile - calculate density profile.
  !----------------------------------------------------------------------  
   	subroutine	density_profile(g_pos, n,bins,min_bin,dens,bin_edges)
   		integer, intent(in)	::n					!Number of objects
   		integer, intent(in)	::bins				!Number of bins to use in histogram
   		real, intent(in)	::g_pos(3,n)		!Positions of objects in group
   		real, intent(in)	::min_bin			!Radius at which to start binning.

   		real, intent(out)	::bin_edges(bins+1)	!Positions of the bin edges
   		real, intent(out)	::dens(bins)		!Density within the bins.
   		
   		real				::r(n)				!Radius of each object from centre
   		real				::vol(bins)			!Volume of each bin
   		integer				::i					!Iterator.

       ! Calculate radii from c_o_m/c_o_d (to make density profiles)
         do i=1,n
           r(i) = sqrt(g_pos(1,i)**2 + g_pos(2,i)**2+ g_pos(3,i)**2)
         end do
         
       ! Calculate Density Profile.
         call LogBins(r,bins,min_bin,dens,bin_edges)
         dens = dens*massarr(2)
         
         do i=1,bins
           vol(i) = (4.0*3.14159*(bin_edges(i+1)**3-bin_edges(i)**3))/3.0
         end do
         
         dens(:) = dens(:)/vol(:)
         if (loud) then
         	 write(*,*)"Calculated spherically symmetric density profile"
         end if
     end subroutine

  !----------------------------------------------------------------------
 	  ! centre_of_mass_offset_fof - calculate the centre of mass
 	  ! offset for a friends-of-friends group
  !----------------------------------------------------------------------  
   	subroutine	centre_of_mass_offset_fof(g_pos, n, mass_offset)
   		integer, intent(in)	::n				!Number of objects in group	
   		real, intent(in)	::g_pos(3,n)	!Co-ords of objects in group
  		
   		real, intent(out)	::mass_offset	!The centre of mass offset
   		
  		real				::most_bound(3)	!The position of the most bound particle
  		real				::c_o_m(3)		!The centre of mass (average position of particles)
  		real				::r_vir			!The virial radius.
  		integer				::i				!Iterator


  		!Find the position of the most bound particle (first particle)
        most_bound = g_pos(:,1)
        
        !Make sure the group has been centred previously.
        !if(most_bound.NE.0.0)then
        !	write(*,*) "ERROR: GROUP WAS NOT CENTRED BEFORE FINDING OFFSETS"
        !	stop
        !end if
        
        	
		!Find the centre of mass of the group.      	
		c_o_m(:) = 0.000
		do i = 1,n
		  c_o_m(1) = c_o_m(1) + g_pos(1,i)
		  c_o_m(2) = c_o_m(2) + g_pos(2,i)
		  c_o_m(3) = c_o_m(3) + g_pos(3,i)
		end do
		c_o_m(:) = c_o_m(:)/real(n)

        
		!Find the virial radius.
		r_vir = (3.d0*n*massarr(2)/(4.d0*pi*delta_vir*rho_c))**(1./3.)
		mass_offset = ((most_bound(1)-c_o_m(1))**2+&
					   &(most_bound(2)-c_o_m(2))**2+&
					   &(most_bound(3)-c_o_m(3))**2)/r_vir
		
     end subroutine
     
  !----------------------------------------------------------------------
 	  ! centre_of_mass_offset_sphere - calculate the centre of mass
 	  ! offset for a spherical overdensity.
  !----------------------------------------------------------------------  
   	subroutine	centre_of_mass_offset_sphere(most_bound, mass_offset)
  		real, intent(in)	::most_bound(3)	!The position of the most bound particle
  		real, intent(out)	::mass_offset	!The centre of mass offset

  		real				::c_o_m(3)		!The centre of mass (average position of particles)
  		real				::r_vir			!The virial radius.
  		integer				::i				!Iterator
  		real				::radius(npart(2))	!The radius from centre of group of all particles in sim.
  		real				::index_array(npart(2)) !The index of all particles 
  		integer				::particles_in_group	!The size of the spherical group.
  		real(8)				::dens			!Density of sphere within current particle

        !ReCentre the whole simulation
        call ReCentre(pos,BoxSize,most_bound)
        
        !Find the radius to ALL PARTICLES from the centre particle of current group
        do i=1,npart(2)
        	radius(i) = sqrt(pos(1,i)**2 + pos(2,i)**2 + pos(3,i)**2)
        end do
        
        !Set an array to index each particle from its original position.
        do i=1,npart(2)
        	index_array(i) = real(i)
        end do
        
        !Sort the radii in ascending order, whilst ordering the indices in the same manner.
        call sort2(npart(2),radius,index_array)
       
        
        !For each particle, estimate the enclosed density, checking whether
        !it is still greater than delta_vir*rho_c.
        do i=2,npart(2)
        	dens = 3*massarr(2)*real(i)/(4*pi*radius(i)**3)
        	!write(*,*) radius(i)
        	if (dens.Lt.delta_vir*rho_c)then
        		particles_in_group = i-1
        		r_vir = radius(i-1)
        		exit
        	end if
        end do
        
		!Find the centre of mass of the group.      	
		c_o_m(:) = 0.000
		do i = 1,particles_in_group
		  c_o_m(1) = c_o_m(1) + pos(1,int(index_array(i)))
		  c_o_m(2) = c_o_m(2) + pos(2,int(index_array(i)))
		  c_o_m(3) = c_o_m(3) + pos(3,int(index_array(i)))
		end do
		c_o_m(:) = c_o_m(:)/real(particles_in_group)
		
		!Calculate the mass_offset
		mass_offset = (c_o_m(1)**2+c_o_m(2)**2+c_o_m(3)**2)/r_vir
		
		!ReCentre the whole dataset back to original co-ords
		call ReCentre(pos,BoxSize,-most_bound)
     end subroutine
     
      
  !----------------------------------------------------------------------
 	  ! read_subgroupv - read a subcat catalogue.
  !---------------------------------------------------------------------- 
  ! npart -- number of particles in (sub)groups
  ! nsub -- number of subhaloes
  ! sublen -- array -- length of each subhalo
  ! suboffset -- array -- offset of each subhalo in IDs catalogue
  ! subparent -- array -- fofcat ID of parent halo

    subroutine read_subgroupv(fname)
      implicit none

      character*(*), intent(in)		:: fname			!input file name
      
      integer 					 	:: i,k				!Iterators
      integer, allocatable			:: sublen(:)       	!Length of each subhalo. Allocated as nsub.
      integer, allocatable			:: suboffset(:)    	!Offset of each subhalo. Allocated as nsub.
      integer, allocatable		 	:: parent_groups(:)	!Temporary array to store which group each subgroup particle belongs to.
      integer					 	:: offset_pos(1)	!Largest subgroup identifier for a particlur subgroup
      integer						:: group_offset		!Offset of the current largest subgroup
      integer						:: group_ID			!ID of current particle belonging to the subgroup.
      
      
      !Open the subgroup file
      open(1,file=fname,status='unknown',form='unformatted')
      
      !Read the number of subgroups.
      read(1) nsub
      
      !Since subparent is a module variable, check if it is already allocated.
      if(allocated(subparent))then
      	  deallocate(subparent)
      end if
      
      !Allocate necessary arrays.
      allocate(sublen(nsub))
      allocate(suboffset(nsub))
      allocate(subparent(nsub))
      
      !Read into the arrays.
      read(1) (sublen(i),i=1,nsub)
      read(1) (suboffset(i),i=1,nsub)
      read(1) (subparent(i),i=1,nsub)
      close(1)

      write(*,*) 'Read ',nsub,' subgroups...'

      allocate(parent_groups(nsub))
      allocate(subgroup_number(npart(2)))
      
      !Initialize Subgroup number to zero, representing belonging to no subgroup.
      subgroup_number = 0
      
      !Now, similarly to read_groups, we cycle through every subgroup, and
      !every particle within that subgroup (from largest to smallest), flagging
      !each particle with a number indicating which subgroup it is in, and also
      !saving a list of numbers referring to which parent group each subgroup is
      !in. Note that the numbers will be fairly randomly dispersed.
      do k = 1,nsub
        offset_pos = maxloc(sublen)
        group_offset = suboffset(offset_pos(1))
      
        parent_groups(k) = subparent(offset_pos(1))
        
        do i=1,maxval(sublen)
          group_ID = ids(group_offset+i)
          
          subgroup_number(group_ID) = k
      
        end do
        
        sublen(offset_pos(1)) = -1
      end do
      
      subparent = parent_groups
    end subroutine

  !----------------------------------------------------------------------
 	  ! number_in_shell - finds how many cartesian positions are in a shell
  !----------------------------------------------------------------------
  	  subroutine number_in_shell(inner_radius, outer_radius, n_objects,g_pos, N)
  	  	real(8), intent(in)		:: inner_radius		!Radius of inner edge of shell
  	  	real(8), intent(in)		:: outer_radius		!Radius of outer edge of shell
  	  	integer, intent(in)		:: n_objects		!Number of objects
  	  	real, intent(in)		:: g_pos(3,n_objects)!Cartesian positions of objects
  	  	
  	  	integer, intent(out)	:: N				!Number of objects in the shell
  	  	
  	  	real(8)					:: ir2, or2, radius
  	  	integer					:: i				!Iterator
  	  	
  	  	ir2 = inner_radius**2
  	  	or2 = outer_radius**2
  	  	N=0
  	  	do i=1,n_objects
  	  		radius = g_pos(1,i)**2 + g_pos(2,i)**2 + g_pos(3,i)**2	
  	  		if(radius.GT.ir2 .AND. radius.LT.or2)then
  	  			N = N+1
  	  		end if
  	  	end do
  	  	
  	  end subroutine
   
  !----------------------------------------------------------------------
 	  ! structure_test - a nonparametric test of structure in shells
  !---------------------------------------------------------------------- 
    !subroutine structure_test(n,g_pos,nside)
    !	integer, intent(in)	:: n			!Number of objects in the group
    !	real, intent(in)	:: g_pos(3,n)	!Positions of objects in group
    !	integer, intent(in)	:: nside		!Resolution parameter for healpix
    !	
    !	integer				:: i			!Iterator
    !	integer				:: pixel(n)		!The pixel occupied by a cartesian vector
    !	real				:: heal_loc(3,n)!The cartesian co-ords of the healpix centres of all group particles	
    !	real				:: radius(n)	!The radius from centre of all particles.
    !	
    !	!What this subroutine does:
    !	!
    !	!I can think of a few ways of doing this. I don't know if they are 
    !	!equivalent or if any is better than any other. The overall idea is to 
    !	!split the group into radial shells, and then for each shell, loop
    !	!through several angular separations, placing down many circles of the
    !	!separation, calculating the density within those volumes, then 
    !	!calculate a variance.
    !		
    !	!I think perhaps the best way to do this is to generate a list of equal
    !	!spaced centres via healpix, and then place my circles around these.
    !		
    !	!Another way would be to use healpix to gather the particles into their
    !	!closest pixel, and do this for several "nside" parameters. May be more
    !	!difficult to calculate volume though.
    !		
    !	!Find the radii of all particles.
    !	do i=1,n
    !		radius(i) = sqrt(g_pos(1,i)**2 + g_pos(2,i)**2 + g_pos(3,i)**2)
    !	end do
    !	
    !	!Snap object positions to healpix centres.
    !	do i=1,n
    !		call vec2pix_nest(nside, g_pos(:,i), pixel(i))
    !	end do
    !	
    !	
    !	!Calculate number of objects assigned to each pixel.
    !	allocate(numbers_of_objects(maxval(pixel)+1)
    !	numbers_of_objects = 0
    !	
    !	do i=1,n
    !		numbers_of_objects(pixel(i)+1) = numbers_of_objects(pixel(i)+1) + 1
    !	end do
    !	
    !	!Convert each number into a density.
    !	vol = some_rubbish
    !	numbers_of_objects = numbers_of_objects*massarr(2)/vol
    !	
    !	
    !			
    !end subroutine
    
  !----------------------------------------------------------------------
 	  ! structure_test_2 - a nonparametric test of structure in shells
  !---------------------------------------------------------------------- 
    subroutine structure_test_2(n,g_pos,n_radial_bins, n_angles,nside,variance)
    	integer, intent(in)	:: n			!Number of objects in the group
    	real, intent(in)	:: g_pos(3,n)	!Positions of objects in group
    	integer, intent(in)	:: n_radial_bins!The number of radial shells to use
    	integer, intent(in)	:: n_angles		!The number of circle sizes to use
    	integer, intent(in)	:: nside		!Resolution parameter for healpix
    	real(8), intent(out):: variance(n_radial_bins,n_angles)	!The density variance
    	
    	integer				:: i,j,k,h		!Iterator
    	real(8), allocatable:: pix_vec(:,:)!The cartesian co-ords of the healpix centres of all group particles	
    	real				:: radius(n)	!The radius from centre of all particles.
    	integer				:: n_pix		!Number of healpix pixels used.
    	real(8)				:: r_vir		!Virial radius of current group
    	real(8)				:: angle, circle_angle, angle_step, step
    	integer	, allocatable:: numbers_of_objects(:)
    	real(8)				:: inner_radius, outer_radius, mean_of_square
    	real(8)				:: square_of_mean, vol, mean
    	real(8), allocatable:: density(:)
    	
    	character(len=255)	:: cwd
    	!What this subroutine does:
    	!
    	!I can think of a few ways of doing this. I don't know if they are 
    	!equivalent or if any is better than any other. The overall idea is to 
    	!split the group into radial shells, and then for each shell, loop
    	!through several angular separations, placing down many circles of the
    	!separation, calculating the density within those volumes, then 
    	!calculate a variance.
    		
    	!I think perhaps the best way to do this is to generate a list of equal
    	!spaced centres via healpix, and then place my circles around these.
    		
    	!Another way would be to use healpix to gather the particles into their
    	!closest pixel, and do this for several "nside" parameters. May be more
    	!difficult to calculate volume though.
    		
    		
    		
    		
    	!Find the radii of all particles.
    	do i=1,n
    		radius(i) = sqrt(g_pos(1,i)**2 + g_pos(2,i)**2 + g_pos(3,i)**2)
    	end do
    	
    	!Output healpix vectors.
    	n_pix = 12*nside**2
    	allocate(pix_vec(3,n_pix))
    	do i=0,n_pix-1
			call pix2vec_nest(nside, i, pix_vec(:,i+1))
			
		end do

    	!Set up values for the loop.
    	r_vir = (3.d0*n*massarr(2)/(4.d0*pi*delta_vir*rho_c))**(1./3.)
    	step = r_vir/n_radial_bins
    	angle_step = (pi/2.)/n_angles
    	variance = 0.d0
    	
    	!Loop over all radial shells out to the virial radius.
    	do k = 1,n_radial_bins
    		inner_radius = (k-1)*step
    		outer_radius = k*step
    		
    		!Loop over angle scales up to 90 degrees.
    		do h = 1, n_angles
    			circle_angle = h*angle_step
    			
				allocate(numbers_of_objects(n_pix))
				numbers_of_objects = 0
				
				do i=1,n_pix
					do j=1,n
						!Find the angle between each pixel and each particle
						angle = acos((pix_vec(1,i)*g_pos(1,j) +&
								& pix_vec(2,i)*g_pos(2,j) +&
								& pix_vec(3,i)*g_pos(3,j))/&
								&sqrt((pix_vec(1,i)**2+pix_vec(2,i)**2+&
								&pix_vec(3,i)**2)*&
								&(g_pos(1,j)**2+g_pos(2,j)**2+g_pos(3,j)**2)))
						if(angle.LT.circle_angle)then
							numbers_of_objects(i) = numbers_of_objects(i) + 1
						end if
					end do
				end do
    	
    			!Convert each number into a density.
    			vol = (2./3.)*pi*(1.d0 - cos(circle_angle))*&
    					&(outer_radius**3 -inner_radius**3)
    					
    			allocate(density(n_pix))
    			density = numbers_of_objects*massarr(2)/vol
    			deallocate(numbers_of_objects)
    			!Find the density variance
    			mean = sum(density)/size(density)
    			!write(*,*) "Min Dens: ", minval(density)
    			!write(*,*) "Max Dens:", maxval(density)
    			!write(*,*) "Mean:", mean
    			do i=1,n_pix
    				variance(k,h) = variance(k,h) + (density(i) - mean)**2
    			end do
    			variance(k,h) = variance(k,h)/n_pix
    			!write(*,*) "Variance:", variance(k,h)
    			!write(*,*) "Densities:", density
    			write(*,*)
    			!square_of_mean = (sum(numbers_of_objects)/&
    			!					&size(numbers_of_objects))**2
    			!mean_of_square = sum(numbers_of_objects**2)/&
    			!					&size(numbers_of_objects)
    			deallocate(density)
    			!variance(k,h) = mean_of_square - square_of_mean
    		end do
    	end do    	
	end subroutine
end module
     

