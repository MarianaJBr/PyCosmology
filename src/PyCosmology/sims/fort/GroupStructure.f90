 
module GroupStructure
	use pix_tools
	use SimOps
	implicit none
	
	
	contains
	
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
    
  !=============================================================================
 	  ! structure_test - a nonparametric test of structure in shells
  !=============================================================================
    subroutine structure_test(n,g_pos,n_radial_bins, n_angles,nside&
    							&,variance,poisson,rads_c,rads,angles_c)
    							
        !What this subroutine does:
    	!
    	!The overall idea is to split the group into radial shells, and then for
    	!each shell, loop through several angular separations, placing down many 
    	!circles of the separation, calculating the density within those 
    	!volumes, then calculate a variance. In effect, this is done in two
    	!ways to keep the volumes equal - by squahsing radial bins as one moves
    	!outwards, and by squashing the angular sizes as one moves outwards.
    		
    	!These are two extremes, and perhaps it would be better to take an 
    	!optimum value within.
    		

     !----------------- Variable Declaration ---------------------------------- 		
    	integer, intent(in)	:: n					!Number of objects in the group
    	real, intent(in)	:: g_pos(3,n)			!Positions of objects in group
    	integer, intent(in)	:: n_radial_bins		!The number of radial shells to use
    	integer, intent(in)	:: n_angles				!The number of circle sizes to use
    	integer, intent(in)	:: nside				!Resolution parameter for healpix
    	
    	real(8), intent(out):: rads_c(n_radial_bins)			!The centres of radial bins with constant step
    	real(8), intent(out):: rads(n_radial_bins)				!The centres of radial bins with diminishin step
    	real(8), intent(out):: angles(n_radial_bins,n_angles)	!The diminishing angles (zeroth angle)
    	real(8), intent(out):: variance(2,n_radial_bins,n_angles)!The density variance
    	real(8), intent(out):: poisson(2,n_radial_bins,n_angles)!The expected variance due to poisson noise for each segment.
    	
    	integer				:: i,j,k,h				!Iterator
    	real(8), allocatable:: pix_vec(:,:)			!The cartesian co-ords of the healpix centres of all group particles	
    	real				:: radius(n)			!The radius from centre of all particles.
    	integer				:: n_pix				!Number of healpix pixels used.
    	real(8)				:: r_vir				!Virial radius of current group
    	real(8)				:: angle				!The angle between an object and healpix pixel.
    	real(8)				:: circle_angle			!The current angle defining a segment
    	real(8)				:: angle_step			!The linear step between angles in the first shell
    	real(8)				:: step					!The linear radial step (also the first radial bin size)
    	integer, allocatable:: numbers_of_objects(:)!The numbers of objects in the current segment size
    	real(8)				:: inner_radius			!The radius defining the inner border of the current bin
    	real(8)				:: outer_radius			!The radius defining the outer border of the current bin
    	real(8)				:: vol					!The volume of the current segment size
    	real(8)				:: mean					!The mean number of particles in the segments
    	real(8), allocatable:: density(:)			!The density of each segment
    	real(8)				:: c_vol_1(n_angles)	!Constant volumes for each angle
    	real(8)				:: V					!Volume of current shell (not just segment)
    	integer				:: mode					!Pancake or pencil-beam
    	real(8)				:: biggest_angle		!The biggest angle allowable in the first shell
    	real(8)				:: smallest_angle		!The smallest angle wanted in the last shell
    	real(8)				:: largest_vol			!The volume that smallest_angle makes in the last shell
    	real(8)				:: ratio				!A number serving a few purposes
    	integer				:: particles			!Number of particles in a shell
    	real(8)				:: mean_particles		!The mean expected particles in each segment.
    	

     !------------------- Variable Initialization ----------------------------
    		
	
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

    	
    	variance = 0.d0
    	poisson = 0.d0
    	
    	
     !-------------------- 3 nested loops ------------------------------------
     		
    	!Use two volume-equalling modes (pencil-beam and pancake)
    	do mode=1,2
    		write(*,*) "mode: ", mode
    		
    		!We need to make sure the radial bins extend to r_vir
    		!This is described by R_1 = r_vir/n^(1/3) - ie. we need 64 bins to reach
    		!down to 0.25 as the first bin. To reach down to 0.1, we need 1000 bins.
    		if(mode==2)then
    			step = r_vir/(real(n_radial_bins)**(1./3.))
    		else
    			step = r_vir/n_radial_bins
    		end if
    		
    		!Set up constant volumes changing radial bins
    		V = (4*pi/3)*step**3
    		
    		!We need to make sure that at least the largest angle in the outmost
    		!radial bin samples the whole space.
    		smallest_angle = acos(pix_vec(1,1)*pix_vec(1,2) + &
    							& pix_vec(2,1)*pix_vec(2,2)&
    							& + pix_vec(3,1)*pix_vec(3,2))/2.d0
    							
    		if(mode==1)then
    			largest_vol = 2*pi*(n_radial_bins**3 -(n_radial_bins-1)**3)*&
    							&(step**3)*(1.d0-cos(smallest_angle))/3.d0   
    			ratio = 3.d0*largest_vol/(2*pi*step**3)
    			if(ratio.GT.2.d0)then
    				write(*,*) "WARNING: You are not sampling the whole space at high radius"
    				biggest_angle = pi/1.2d0
    			else
    				biggest_angle = acos(1.d0-3.d0*largest_vol/(2*pi*step**3))
    			end if
    			angle_step = (biggest_angle)/n_angles
    		else
    			angle_step = 1.5d0*smallest_angle/n_angles
    		end if
        	outer_radius = 0.d0
    		
    		!Set up constant volumes changing angles
    		if(mode==1)then
    			do i=1,n_angles
    				circle_angle =i*angle_step
    				c_vol_1(i) = (2./3.)*pi*(1.d0 - cos(circle_angle))*step**3
    			end do
    		end if
    		
    		!Loop over all radial shells out to the virial radius.
    		do k = 1,n_radial_bins
    			if(mode==1)then
    				inner_radius = outer_radius
    				outer_radius = inner_radius + step
    				rads_c(k) =((0.5*(inner_radius**3+outer_radius**3))**(1./3.))/&
    							&r_vir
    			elseif(mode==2)then
    				if(k==1)then
    					outer_radius = 0.d0
    				end if
    				inner_radius = outer_radius
    				outer_radius = (3*V/(4*pi) + inner_radius**3)**(1./3.)

    				rads(k) =((0.5*(inner_radius**3 + outer_radius**3))**(1./3.))/&
    							&r_vir

    			end if
    			
    			!Loop over angle scales up to 90 degrees.
    			do h = 1, n_angles
    				if(mode==1)then
    					!Modify angle to keep constant volume
    					ratio = 1.d0-3.d0*c_vol_1(h)/(2*pi*&
    											&(outer_radius**3-inner_radius**3))
    					if(ratio.LT.-1.d0)then
    						ratio = -1.d0
    					end if
    					circle_angle = acos(ratio)

    				else
    					!Angle increases linearly
    					circle_angle = h*angle_step
    					angles_c(h) = circle_angle
    				end if
    				
					allocate(numbers_of_objects(n_pix))
					numbers_of_objects = 0
					
					vol = (2./3.)*pi*(1.d0 - cos(circle_angle))*&
    						&(outer_radius**3 -inner_radius**3)
    						
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
								if(radius(j).GT.inner_radius.AND.radius(j).LT.outer_radius)then
									numbers_of_objects(i) = numbers_of_objects(i) + 1
								end if
							end if
						end do
					end do
    		
					!Find the Poisson noise in the current shell
						!First find the number of objects in shell
					particles = 0
					do i=1,n
						if(radius(i).GT.inner_radius.AND.radius(i).LT.outer_radius)then
							particles = particles + 1
						end if
					end do
					!Calculate ratio of volume taken by segment
					ratio = circle_angle/(2*pi)
					
					!Find mean number in each segment
					mean_particles = ratio*particles
					
					!Then variance = mean for Poisson.
						
					poisson(mode,k,h) = mean_particles*massarr(2)/vol
					
    				!Convert each number into a density.
    				allocate(density(n_pix))
    				density = numbers_of_objects*massarr(2)/vol
    				deallocate(numbers_of_objects)
    				
    				!Find the density variance
    				mean = sum(density)/size(density)
    				do i=1,n_pix
    					variance(mode,k,h) = variance(mode,k,h)+(density(i)-mean)**2
    				end do
    				variance(mode,k,h) = variance(mode,k,h)/n_pix
        	
    				deallocate(density)
        	
    			end do
    		end do  
    	end do
	end subroutine
end module
