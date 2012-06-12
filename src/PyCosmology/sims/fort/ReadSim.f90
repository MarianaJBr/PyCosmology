
module SimOps
	implicit none


	real, allocatable		:: pos(:,:),vel(:,:)
	real, allocatable		:: g_pos(:,:), g_vel(:,:)
	integer,allocatable  	:: ID(:)
	integer,allocatable  	:: grouplen(:), groupoffset(:), ids(:),group_IDS(:)
	integer, allocatable	:: group_number(:)
	
	real, allocatable 		:: pos_r(:,:), r(:), vel_r(:,:)
	real(8),allocatable   	:: dens(:), bin_edges(:), vol(:)
	
    
	!Variables defined in the header of a snapshot file.	
	integer  :: FlagSfr, FlagFeedback, Nall(6), FlagAge, FlagCooling 
	integer  :: FlagMetals, NallHW(6), FlagEntrICs, npart(6),NumFiles

	real(8)  ::massarr(6), time, BoxSize,OmegaLambda
	real(8)	 ::HubbleParam,Omega0, Redshift
	
	integer		:: ngroups
	real		:: axis_ratio_b, axis_ratio_c
	logical		:: loud = .FALSE.
	
	contains


  !----------------------------------------------------------------------
 	  ! ReadSim - Simply read in information from the simulation
  !----------------------------------------------------------------------		
    subroutine ReadSim(filename)    
		integer :: i
		character(len=*), intent(in) :: filename

		! Read in Simulation File values (only positions, velocities and ID's)  
		open(unit=10,file=filename,form='UNFORMATTED')
        
		read(10) npart, massarr, time, Redshift, FlagSfr, FlagFeedback, Nall,&
				  & FlagCooling, NumFiles, BoxSize, Omega0, OmegaLambda,&
				  &HubbleParam,FlagAge, FlagMetals, NallHW, FlagEntrICs
																							
		write(*,*) "Number of particles(1) = ", npart(1)
		write(*,*) "Number of particles(2) = ", npart(2)
		write(*,*) "Number of particles(3) = ", npart(3)
		write(*,*) "Number of particles(4) = ", npart(4)
		write(*,*) "Number of particles(5) = ", npart(5)
		write(*,*) "Number of particles(6) = ", npart(6)
		write(*,*) "Time of snapshot: ", time
		write(*,*) "Redshift of snapshot: ", Redshift
		write(*,*) "BoxSize of Sim: ", BoxSize
		
		allocate(pos(3,npart(2)),vel(3,npart(2)),ID(npart(2)))
		read(10)
		read(10)
				
		
		read(10) (ID(i), i=1,npart(2))
		write(*,*) "Read Particle ID's"
		rewind(10)
		read(10)
		read(10) (pos(1,ID(i)),pos(2,ID(i)),pos(3,ID(i)), i=1,npart(2))
		write(*,*) "Read Positions"
		read(10) (vel(1,ID(i)),vel(2,ID(i)),vel(3,ID(i)), i=1,npart(2))
		write(*,*) "Read Velocities"
		close(10)
      
      deallocate(ID)
    end subroutine
    
  !----------------------------------------------------------------------
 	  ! find_groups - establish which objects belong to which groups.
  !----------------------------------------------------------------------
    subroutine find_groups(filename,ids_file)  
    
    	integer 			:: k, ntot, group_offset, group_size
    	character(len=*)	:: filename,ids_file
    	
    	integer				::offset_pos(1), group_ID, i

    	
    	
    	open(1,file=filename,status='unknown',form='unformatted')
    	          
    	read(1) ngroups                         ! Total number of FOF groups
    	
    	write(*,'(a,i6,a)') 'Reading ',ngroups,' groups...'
    	
    	allocate(grouplen(ngroups),groupoffset(ngroups))
    	
    	read(1) (grouplen(i),i=1,ngroups)      ! Length of each FOF group
    	write(*,*) 'Read group lengths'
    	read(1) (groupoffset(i),i=1,ngroups)  ! Offset in the pos.ids file
    	write(*,*) 'Read group offsets'
    	!read(1) (nsubgroups(i),i=1,ngroups) ! Number of subgroups -- ignore
    	!write(*,*) 'Read number of subgroups'
    	close(1)
    	
    	open(1,file=ids_file,status='unknown',form='unformatted')
    	read(1) ntot                     ! Total number of particles in FOF groups
    	allocate(ids(ntot))
    	read(1) (ids(i),i=1,ntot)  ! Particle IDs
    	close(1)    

        allocate(group_number(size(pos(1,:))))
        group_number(:) = 0
    
    ! Locate Largest Group/Halo and save properties.
        do k = 1,ngroups
        	if(loud)then
        		write(*,*) "Calculating characteristics of group", k
        	end if
          offset_pos = maxloc(grouplen)
          group_offset = groupoffset(offset_pos(1))
          group_size = maxval(grouplen)
          if (loud)then
          	  write(*,*) "Group Size: ", group_size
          end if



          if (loud) then
          	  write(*,*) "Specifying group particles"
          end if
          
          do i=1,maxval(grouplen)
            group_ID = ids(group_offset+i)
            
            group_number(group_ID) = k
     
          end do
          
          grouplen(offset_pos(1)) = -1
        end do
      end subroutine
 
  !----------------------------------------------------------------------
 	  ! centre - change co-ords of a group to its centre
  !----------------------------------------------------------------------        
        
  	subroutine centre(g_pos,n,mode)
  	!
  	!mode refers to where to place the centre (0,0,0).
  	!mode = 1 uses the first particle (centre of density)
  	!mode = 2 uses the average particle position (centre of mass)
  	!mode = 3 recentres with equal limits in both directions.
  	!
  	
    	integer,intent(in)	:: n 
    	integer,intent(in),optional :: mode
  		real,intent(inout)	:: g_pos(3,n)
  		
  		real				:: centre_of_density(3), c_o_m(3), equal_centre(3)
  		integer				:: i, temp_mode
 		
  		if (present(mode)) then
  			temp_mode = mode
  		else
  			temp_mode = 1
  		end if
  		

  		
        centre_of_density = g_pos(:,1)
        ! Re-centre the group around the centre of density
        call ReCentre(g_pos, real(BoxSize),centre_of_density)
        
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
				call ReCentre(g_pos, real(BoxSize),c_o_m)
		end if
        
		if (temp_mode == 3) then
			equal_centre(1) = (maxval(g_pos(1,:))+minval(g_pos(1,:)))/2.d0
			equal_centre(2) = (maxval(g_pos(2,:))+minval(g_pos(2,:)))/2.d0
			equal_centre(3) = (maxval(g_pos(3,:))+minval(g_pos(3,:)))/2.d0
			! Re-centre the group around the centre
				call ReCentre(g_pos, real(BoxSize),equal_centre)
		end if
		
		if (loud) then
        write(*,*) "Group Re-centred"
        end if
     end subroutine
     
            
  
  !----------------------------------------------------------------------
 	  ! rotate - rotate data into moment of inertia defined coords plus angle
  !----------------------------------------------------------------------            
	subroutine rotate(g_pos, g_vel, n, u, theta)
		integer, intent(in)	::n
		real, intent(inout)	::g_pos(3,n), g_vel(3,n)
		real, intent(in)	::u(3),theta

		
		integer	   :: nrot, i
    	real      :: IT(3,3), eigval(3), pos_r(3,n)
		real      :: eigvec(3,3),bin_size

		
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
              if (loud) then
              write(*,*) "Ratio of axes lengths b/a = ", axis_ratio_b
              write(*,*) "Ratio of axes lengths c/a = ", axis_ratio_c
              end if
              
            ! Project data into co-ordinates defined by
                ! xbar = eigvec(1)*(x,y,z)
                ! ybar = eigvec(2)*(x,y,z)
                ! zbar = eigvec(3)*(x,y,z)
                  
              do i=1,n
                g_pos(:,i) = matmul(eigvec,g_pos(:,i))
                g_vel(:,i) = matmul(eigvec,g_vel(:,i))
              end do
              
              call GenRotate(g_pos,u,theta,pos_r)
              g_pos = pos_r

              call GenRotate(g_vel,u,theta,pos_r)
              g_vel = pos_r
      end subroutine
      
      
  !----------------------------------------------------------------------
 	  ! density_profile - calculate density profile.
  !----------------------------------------------------------------------  
   	subroutine	density_profile(g_pos, n,bins,min_bin,dens,bin_edges)
   		integer, intent(in)	::n, bins
   		real, intent(in)	::g_pos(3,n),min_bin

   		real, intent(out)				::bin_edges(bins+1),dens(bins)
   		real							::r(n), vol(bins)
   		integer							::i

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
         
  
  subroutine ReCentre(pos, BoxSize, centre)
    real, intent(in)    :: BoxSize
    real, intent(inout) :: pos(:,:)
    real,intent(in)      :: centre(3)
    integer             :: i,j
    
    
    do j=1,size(pos(1,:))
      do i=1,size(pos(:,1))
        if (pos(i,j).LT.centre(i)-BoxSize/2.0)then
          pos(i,j) = pos(i,j) + BoxSize
        else if(pos(i,j).GT.centre(i)+BoxSize/2.0)then
          pos(i,j) = pos(i,j) - BoxSize
        end if
        
        
      end do
      pos(:,j) = pos(:,j)-centre(:)
    end do
  end subroutine
  
  
  function Inertia(mass,pos)
    real, intent(in)  :: mass, pos(:,:)
    real, allocatable :: pos_1(:,:)
    integer           :: i,j,k, length
    real, allocatable :: r(:)
    real              :: c_o_m(3)
    real              :: Inertia(3,3)
    
    length = size(pos(1,:))
    c_o_m(:) = 0.000
    do i = 1,length
      c_o_m(:) = c_o_m(:) + pos(:,i)
    end do
    c_o_m = c_o_m/real(length)
    
    allocate(r(length))
    allocate(pos_1(3,length))
    do i = 1,length
      pos_1(:,i) = pos(:,i) - c_o_m(:)
    end do
    
    do i =1,length
      r(i) = sqrt(pos_1(1,i)**2 + pos_1(2,i)**2 + pos_1(3,i)**2)
    end do
    
    
    Inertia(:,:) = 0.0
    do j=1,3
      do k=1,3
        do i=1,length
          if (j==k) then
            Inertia(j,k) = Inertia(j,k) + (r(i)**2- pos_1(j,i)*pos_1(k,i))
          else
            Inertia(j,k) = Inertia(j,k) - pos_1(j,i)*pos_1(k,i)
          end if
        end do
      end do
    end do
    Inertia = mass*Inertia

  end function
  
  
  subroutine EigOrder(eigval,eigvec)
    real, intent(inout) :: eigval(:), eigvec(:,:)
    real                :: t_eigval(size(eigval))
    real                :: t_eigvec(size(eigvec(:,1)),size(eigvec(1,:)))
    integer             :: maxpos(1)
    integer :: i
    
    do i=1,3
      maxpos = maxloc(eigval)
      t_eigval(i) = eigval(maxpos(1))
      eigval(maxpos(1)) = minval(eigval)-100.00
      
      t_eigvec(:,i) = eigvec(:,maxpos(1))
    end do
    
    eigval = t_eigval
    eigvec = t_eigvec
  end subroutine
  
  subroutine GenRotate(pos,u,theta,pos_r)
    real, intent(in)    ::u(:), theta,pos(:,:)
    real, intent(out) :: pos_r(:,:)
    real              :: ROT(3,3)
    integer           :: i
    
    ROT(1,1) = cos(theta) + (u(1)**2)*(1-cos(theta))
    ROT(1,2) = u(1)*u(2)*(1-cos(theta))-u(3)*sin(theta)
    ROT(1,3) = u(1)*u(2)*(1-cos(theta))+u(2)*sin(theta)
    
    ROT(2,1) = u(2)*u(1)*(1-cos(theta))+u(3)*sin(theta)
    ROT(2,2) = cos(theta) + u(2)*u(2)*(1-cos(theta))
    ROT(2,3) = u(2)*u(3)*(1-cos(theta))-u(1)*sin(theta)
    
    ROT(3,1) = u(3)*u(1)*(1-cos(theta))-u(2)*sin(theta)
    ROT(3,2) = u(3)*u(2)*(1-cos(theta))+u(1)*sin(theta)
    ROT(3,3) = cos(theta)+u(3)*u(3)*(1-cos(theta))
    
    do i = 1,size(pos(1,:))
      pos_r(:,i) = matmul(ROT,pos(:,i))
    end do
  end subroutine
    
  subroutine LogBins(dat, bins, min_bin, binned,bin_edges)
    real, intent(in)  :: dat(:), min_bin
    integer, intent(in):: bins
    real, intent(out)   :: binned(bins)
    real                :: bin_edges(:), bin_size
    integer :: i,j
    
    bin_size = (maxval(dat)/min_bin)**(1.0/(bins-1))

    bin_edges(1) = 0.0d0
    do i=1,bins
      bin_edges(i+1) = min_bin*bin_size**(i-1)
    end do
    

    binned(:) = 0.0d0
    do i=1,size(dat)
      do j=1,bins
        if (dat(i).GT.bin_edges(j).AND.dat(i).LT.bin_edges(j+1))then
          binned(j) = binned(j) + 1.0d0
        end if
      end do
    end do
  end subroutine
    
    


    
    
        SUBROUTINE jacobi(a,n,np,d,v,nrot)
      INTEGER n,np,nrot,NMAX
      REAL a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))&
            &*g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      return
      END
  
  
      
    
    


    
end module
     

