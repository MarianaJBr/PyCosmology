module shapelets
  implicit none
  
  !A module defining various subroutines and functions use in 3D shapelet
  !decomposition.
  	  
  !There are several subroutines defining significant components of the
  !process, and one subroutine to join them together. All interesting quantities
  !are stored as module-level variables for access from Python.
    
  !----- Universal Parameters --------------------------------------------------
  	real(16), parameter	   ::pi= 3.141592653589793  !Value of Pi

  !----- Accessible Quantities -------------------------------------------------
 	real(8)					:: beta_c			!Beta parameter from Fluke.
 	real(8)					:: dx				!The interval between cells.
 	integer					:: nmax				!The maximum shapelet order
 	
 	real(8), allocatable	:: Ints(:,:)		!Integral factors. Takes shape (Ng, nmax+1)
 	real(8), allocatable	:: coeffs(:,:,:)	!The cube of coefficients of shapelets.
 	real(8), allocatable	:: recon(:,:,:)		!The reconstructed group density.
 	
 	
 	
 	
 	
 	
 	
 	contains
  !----------------------------------------------------------------------
 	  ! shapelet_driver - drives the other subroutines.
  !---------------------------------------------------------------------- 
  	  subroutine shapelet_driver(density, x_max, do_reconstruct, Ng)
  	  	real(8), intent(in)		:: density(Ng,Ng,Ng)
  	  	real(8), intent(in)		:: x_max
  	  	integer, intent(in)		:: Ng
  	  	logical, intent(in)		:: do_reconstruct
  	  	
  	  
 	  	dx = 2.d0*x_max/Ng
 	  	nmax = (Ng-3)/2
 	  	beta_c = 0.8*x_max/sqrt(2.0d0*Ng)
 	  	
 	  	
 	  	write(*,*) "x_max = ", x_max
 	  	write(*,*) "Ng = ", Ng
 	  	write(*,*) "nmax = ", nmax
 	  	write(*,*) "dx = ", dx
 	  	write(*,*) "beta = ", beta_c
 	    
 	  	
  	  	call coeff_cube(density, Ng, x_max)
 	
  	  	if(do_reconstruct)then
  	  		call reconstruct_cube(Ng, x_max)
  	  		write(*,*) 'f(26,27,27):', density(26,27,27)
  	  		write(*,*) 'recon(26,27,27):', recon(26,27,27)
  	  	end if
 	
  	  end subroutine
  	  
  	  
  !----------------------------------------------------------------------
 	  ! BASIS_VECTOR - calculates bases recursively
  !----------------------------------------------------------------------  
  	  subroutine basis_vector(x,n,bases)
  	  !Calculates the first n+1 basis functions at a point x, based on the
  	  !iterative method of Fluke et. al.
  	  	  
  	  	real(16), intent(in)	::x
  	  	integer, intent(in)		::n
  	  	
  	  	real(8), intent(out)	::bases(n+1)
  	  	integer					::i
  	  	
  	  	bases(1) = 1.d0/(sqrt(beta_c*sqrt(pi)))*exp(-x ** 2.d0/(2.d0*beta_c**2))
  	  	if(n.GT.0)then
  	  		bases(2) = sqrt(2.d0)*x*bases(1)/beta_c
  	  	end if
  	  	if (n.GT.1)then
			do i=3,n+1
				bases(i) = (x/beta_c)*sqrt(2.d0/(i-1))*bases(i-1)-&
				&sqrt(real(i-2)/real(i-1))*bases(i-2)
			end do
		end if
  	  	

  	  end subroutine
  	  

  !----------------------------------------------------------------------
 	  ! cube_ints - finds the integrals necessary for a cube
  !----------------------------------------------------------------------	  
 	subroutine cube_ints(x_max,Ng)
 	!Calculates the first n+1 Integral factors for all Ng cells. Beta must be
 	!specified beforehand.
 	
 	  integer, intent(in)  ::Ng
 	  real(8), intent(in)	::x_max
 	  
 	  integer              ::i,j
 	  real(8)              ::bases(Ng+1,nmax+1)
 	  real(16)				:: x_a, x_b, erf_a, erf_b, a, b

 	  do i=1,Ng+1
 	  	a = -x_max + (i-1)*2.d0*x_max/Ng
 	  	call basis_vector(a,nmax,bases(i,:))
 	  end do
 	  
 	  if(allocated(Ints))then
 	  	  deallocate(Ints)
 	  end if
 	  
 	  allocate(Ints(Ng,nmax+1))
 	  
 	  Ints = 0.0d0

 	  do i=1,Ng
 	  	a = -x_max + (i-1)*2.d0*x_max/Ng
 	  	b = a + 2.d0*x_max/Ng 	 
 	  	x_a = a/(beta_c*sqrt(2.d0))
 	  	x_b = b/(beta_c*sqrt(2.d0))
 	  	
 	  	write(*,*) "a = ", a
 	  	write(*,*) "b = ", b
 	  	write(*,*) "x_a = ", x_a
 	  	write(*,*) "x_b = ", x_b
 	  	
 	  	erf_a = erf(x_a)
 	  	erf_b = erf(x_b)
		write(*,'(A9,es28.20)') "erf(a) = ", log(abs(erf_a))
		write(*,'(A9,es28.20)') "erf(b) = ", log(abs(erf_b))
		
 	  	Ints(i,1) = (beta_c*pi**0.5d0 /2.d0)**(0.5d0) * (erf_b-erf_a)
		  
		Ints(i,2) = -sqrt(2.d0)*beta_c*(bases(i+1,1)-bases(i,1))
 	  
		do j=3,nmax+1
 	    	Ints(i,j) = -beta_c*sqrt(2.d0/(j-1))*(bases(i+1,j-1)-bases(i,j-1))+&
 	                &sqrt((i-1)/real(i))*Ints(i,j-2)
 	    end do
 	  end do
 	  
 	  do i=1,nmax+1
 	  	write(*,'(es28.20)')Ints(1,i)
 	  end do
 	  	
 	end subroutine

  !----------------------------------------------------------------------
 	  ! COEFF_CUBE evaluates the coefficients in the case of a nice cube.
 	  	  ! This follows strictly the algorithm and optimization layed out
 	  	  ! in Fluke et. al.
  !----------------------------------------------------------------------	 	 
 	subroutine coeff_cube(f, Ng,x_max_8)
 	  integer, intent(in)   :: Ng
 	  real(8), intent(in)   :: f(Ng,Ng,Ng)
 	  real(8), intent(in)   :: x_max_8

 	  
 	  integer               :: i,j,k,ii,jj,kk
 	  real(8)				 :: x_max
 	  

 	  if(allocated(coeffs))then
 	  	  deallocate(coeffs)
 	  end if
 	  
 	  allocate(coeffs(nmax+1,nmax+1,nmax+1))
 
 	  coeffs = 0.0d0
 	  call cube_ints(x_max_8,Ng)
 

 	  do k=1,nmax+1
 	  	i=1
 	  	j=1
 	  	if (i+j+k .GT. nmax+1)then
 	      		exit
 	    end if
 	    do j=1,nmax+1
 	    	i=1
 	    	if (i+j+k .GT. nmax+1)then
 	       		exit
 	      	end if
 	      do i=1,nmax+1
 	      	if (i+j+k .GT. nmax+1)then
 	      		exit
 	      	end if
 	      	do kk=1,Ng
 	      	  do jj = 1,Ng
 	      	    do ii = 1,Ng
 	      	    	coeffs(i,j,k) = coeffs(i,j,k) + &
 	      	    	& f(ii,jj,kk)*Ints(ii,i)*Ints(jj,j)*Ints(kk,k)
 	      	    end do
 	      	  end do
 	      	end do
 	      end do
 	    end do
 	  end do
 	end subroutine

  !----------------------------------------------------------------------
 	  ! RECONSTRUCT_CUBE reconstructs the data in the case of a nice cube.
 	  	  ! This follows strictly the algorithm and optimization layed out
 	  	  ! in Fluke et. al.
  !----------------------------------------------------------------------	
  	  subroutine reconstruct_cube(Ng,x_max)
  	  	integer, intent(in)	:: Ng
  	  	real(8), intent(in)	:: x_max

  	  	integer	:: i,j,k,ii,jj,kk
  	  	real(8)	:: bases(Ng,nmax+1)
  	  	real(16)::x
  	  	
  	  	if(allocated(recon))then
  	  		deallocate(recon)
  	  	end if
  	  	
  	  	allocate(recon(Ng,Ng,Ng))
  	  	
  	  	recon = 0.0d0
  	  	
  	  	do j = 1,Ng
  	  		x = -x_max+(j-0.5d0)*dx
  	  		call basis_vector(x,nmax,bases(j,:))
  	  	end do
  	  	
  	  	do k=1,Ng
  	  		do j = 1,Ng
  	  			do i = 1,Ng
  	  				do kk = 1,nmax+1
  	  					jj=1
  	  					ii=1
  	  					!if (ii+jj+kk .GT. nmax+1)then
  	  					!	exit
  	  					!end if
  	  					do jj = 1,nmax+1
  	  						ii=1
  	  					!	if (ii+jj+kk .GT. nmax+1)then
  	  					!		exit
  	  					!	end if	
  	  						do ii = 1,nmax+1
  	  					!	if (ii+jj+kk .GT. nmax+1)then
  	  					!		exit
  	  					!	end if
  	  							recon(i,j,k) = recon(i,j,k) + &
  	  								& coeffs(ii,jj,kk)*bases(i,ii)*&
  	  								&bases(j,jj)*bases(k,kk)
  	  						end do
  	  					end do
  	  				end do
  	  			end do
  	  		end do
  	  	end do
  	  end subroutine
  	  
  	  
  !----------------------------------------------------------------------
 	  ! TEST BASIS_VECTOR - create reconstructions from singular basis
 	  !							functions.
  !----------------------------------------------------------------------  
  	  subroutine test(Ng,cube_0,cube_1,cube_2,cube_3)
  	  	!This subroutine produces four cubes, based on the basis vectors
  	  	!with beta = 1, n=(0,0,0), (0,1,0),(2,0,2),(1,2,4).
  	  		
  	  	integer, intent(in)	::Ng
  	  	real(8), intent(out)::cube_0(Ng,Ng,Ng), cube_1(Ng,Ng,Ng)
  	  	real(8), intent(out)::cube_2(Ng,Ng,Ng), cube_3(Ng,Ng,Ng)
  	  	
  	  	real(8)		:: x_max
  	  	real(8)		:: base(5,Ng)
  	  	real(16)	:: x
  	  	integer		::i,j,k
  	  	
  	  	
  	  	x_max = 1.5d0
  	  	beta_c = 1.d0
  	  	
  	  	cube_0 = 0.0d0
  	  	cube_1 = 0.0d0
  	  	cube_2 = 0.0d0
  	  	cube_3 = 0.0d0
  	  	
  	  	do i=1,Ng
  	  		x = -x_max+(i-0.5d0)*2.0*x_max/Ng
  	  		call basis_vector(x,4,base(:,i))
  	  	end do
  	  	
  	  	do k=1,Ng
  	  		do j = 1,Ng
  	  			do i = 1,Ng
  	  				cube_0(i,j,k) =cube_0(i,j,k) + base(1,i)*base(1,j)*base(1,k)
  	  				cube_1(i,j,k)=cube_1(i,j,k)+base(1,i)*base(2,j)*base(1,k)
  	  				cube_2(i,j,k)=cube_2(i,j,k)+base(3,i)*base(1,j)*base(3,k)
  	  				cube_3(i,j,k)=cube_3(i,j,k)+base(2,i)*base(3,j)*base(5,k)
  	  			end do
  	  		end do
  	  	end do

  	  end subroutine
end module
 	  
