module shapelets
	implicit none
	
  !A module defining various subroutines and functions use in 3D shapelet
  !decomposition.
  	  
  !There are several subroutines defining significant components of the
  !process, and one subroutine to join them together.
  	  
  !----- Universal Parameters --------------------------------------------------
  	real(16), parameter	   ::pi= 3.141592653589793  !Value of Pi

  !----- Accessible Quantities -------------------------------------------------
 	real(8)					:: beta				!Beta parameter from Fluke.
 	real(8)					:: dx				!The interval between cells.
 	integer					:: nmax				!The maximum shapelet order
 	
	contains

  !----------------------------------------------------------------------
 	  ! shapelet_driver - drives the other subroutines.
  !---------------------------------------------------------------------- 
  	  subroutine shapelet_driver(density,x_max,do_reconstruct, Ng,coeffs,recon)
  	  	real(8), intent(in)		:: density(Ng,Ng,Ng)
  	  	real(8), intent(in)		:: x_max
  	  	integer, intent(in)		:: Ng
  	  	logical, intent(in)		:: do_reconstruct
  	  	
  	  	real(8), intent(out)	:: coeffs((Ng-3)/2+1,(Ng-3)/2+1,(Ng-3)/2+1)
  	  	real(8), intent(out)	:: recon(Ng,Ng,Ng)
  	  	
 	  	dx = 2.d0*x_max/Ng
 	  	nmax = (Ng-3)/2
 	  	beta = 0.8*x_max/sqrt(2.0d0*Ng)
 	  	
 	  	
 	  	write(*,*) "x_max = ", x_max
 	  	write(*,*) "Ng = ", Ng
 	  	write(*,*) "nmax = ", nmax
 	  	write(*,*) "dx = ", dx
 	  	write(*,*) "beta = ", beta
 	    
 	  	
  	  	call coeff_cube(density, Ng, x_max,coeffs)
 	
  	  	if(do_reconstruct)then
  	  		call reconstruct_cube(coeffs,Ng, x_max,recon)
  	  		write(*,*) 'f(26,27,27):', density(26,27,27)
  	  		write(*,*) 'recon(26,27,27):', recon(26,27,27)
  	  	end if
 	
  	  end subroutine 	
  	  
  !----------------------------------------------------------------------
 	  ! cyclic_test - performs a basic test of the de/re-construction
  !---------------------------------------------------------------------- 
  	  subroutine cyclic_test(bases,n_b,x_max,Ng,intermediate_recon,coeffs,recon)
  	  	integer, intent(in)		:: n_b
  	  	integer, intent(in)		:: bases(n_b,3)
  	  	real(8), intent(in)		:: x_max
  	  	integer, intent(in)		:: Ng
  	  	
  	  	real(8), intent(out)	:: intermediate_recon(Ng,Ng,Ng)
  	  	real(8), intent(out)	:: coeffs((Ng-3)/2+1,(Ng-3)/2+1,(Ng-3)/2+1)
  	  	real(8), intent(out)	:: recon(Ng,Ng,Ng)

  	  	integer					:: i,j,k
  	  	
 	  	dx = 2.d0*x_max/Ng
 	  	nmax = (Ng-3)/2
 	  	beta = 0.8*x_max/sqrt(2.0d0*Ng)


 	  	coeffs(:,:,:) = 0.d0
 	    do i=1,n_b
 	    write(*,*) i, 'done'
 	    	coeffs(bases(i,1)+1,bases(i,2)+1,bases(i,3)+1) = 1.d0
 	    end do
 	    
 	    
 	  	write(*,*) "x_max = ", x_max
 	  	write(*,*) "Ng = ", Ng
 	  	write(*,*) "nmax = ", nmax
 	  	write(*,*) "dx = ", dx
 	  	write(*,*) "beta = ", beta
 	    
 	  	
  	  	call reconstruct_cube(coeffs,Ng, x_max,recon)
  	  	
  	  	intermediate_recon = recon
  	  	
  	  	call coeff_cube(recon, Ng, x_max,coeffs)
 
  	  	call reconstruct_cube(coeffs,Ng, x_max,recon)
  	  	!Check how the coefficients match up.
 	    do i=1,n_b
 	    	write(*,*) bases(:,i),coeffs(bases(i,1)+1,bases(i,2)+1,bases(i,3)+1)
 	    end do
 	    
 	    write(*,*) "THE FOLLOWING COEFFICIENTS ARE LARGER THAN 0:"
 	    do i=1,Ng
 	    	do j=1,Ng
 	    		do k=1,Ng
 	    			if (coeffs(i,j,k).GT.1E-16) then
 	    				write(*,*) "(",i-1,',',j-1,',',k-1,') :', coeffs(i,j,k)
 	    			end if
 	    		end do
 	    	end do
 	    end do
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
  	  	
  	  	bases(1) = 1.d0/(sqrt(beta*sqrt(pi)))*exp(-x ** 2.d0/(2.d0*beta**2))
  	  	if(n.GT.0)then
  	  		bases(2) = sqrt(2.d0)*x*bases(1)/beta
  	  	end if
  	  	if (n.GT.1)then
			do i=3,n+1
				bases(i) = (x/beta)*sqrt(2.d0/(i-1))*bases(i-1)-&
				&sqrt(real(i-2)/real(i-1))*bases(i-2)
			end do
		end if
  	  	

  	  end subroutine
  	  
  	  
  !----------------------------------------------------------------------
 	  ! cube_ints - finds the integrals necessary for a cube
  !----------------------------------------------------------------------	  
 	subroutine cube_ints(x_max,Ng,Ints)
 	!Calculates the first n+1 Integral factors for all Ng cells. Beta and nmax
 	!must be specified beforehand.
 	
 	  integer, intent(in)  ::Ng
 	  real(8), intent(in)	::x_max
 	  
 	  real(8), intent(out)	:: Ints(Ng,nmax+1)
 	  
 	  integer              ::i,j
 	  real(8)              ::bases(Ng+1,nmax+1)
 	  real(16)				:: x_a, x_b, erf_a, erf_b, a, b

 	  do i=1,Ng+1
 	  	a = -x_max + (i-1)*2.d0*x_max/Ng
 	  	call basis_vector(a,nmax,bases(i,:))
 	  end do

 	  Ints = 0.0d0

 	  do i=1,Ng
 	  	a = -x_max + (i-1)*2.d0*x_max/Ng
 	  	b = a + 2.d0*x_max/Ng 	 
 	  	x_a = a/(beta*sqrt(2.d0))
 	  	x_b = b/(beta*sqrt(2.d0))

 	  	erf_a = erf(x_a)
 	  	erf_b = erf(x_b)

		
 	  	Ints(i,1) = (beta*pi**0.5d0 /2.d0)**(0.5d0) * (erf_b-erf_a)
		  
		Ints(i,2) = -sqrt(2.d0)*beta*(bases(i+1,1)-bases(i,1))
 	  
		do j=3,nmax+1
 	    	Ints(i,j) = -beta*sqrt(2.d0/(j-1))*(bases(i+1,j-1)-bases(i,j-1))+&
 	                &sqrt((i-1)/real(i))*Ints(i,j-2)
 	    end do
 	  end do
 	  	
 	end subroutine
 	
 	
  !----------------------------------------------------------------------
 	  ! COEFF_CUBE evaluates the coefficients in the case of a nice cube.
 	  	  ! This follows strictly the algorithm and optimization layed out
 	  	  ! in Fluke et. al.
  !----------------------------------------------------------------------	 	 
 	subroutine coeff_cube(f, Ng,x_max_8,coeffs)
 	  integer, intent(in)   :: Ng
 	  real(8), intent(in)   :: f(Ng,Ng,Ng)
 	  real(8), intent(in)   :: x_max_8

 	  
 	  integer               :: i,j,k,ii,jj,kk
 	  
 	  real(8)				 :: Ints(Ng,nmax+1)
 	  real(8), intent(out)  :: coeffs(nmax+1,nmax+1,nmax+1)
 	  
 	  call cube_ints(x_max_8,Ng,Ints)

 	  write(*,*) "NOW NMAX IS: ", nmax
 	  coeffs = 0.0d0
 	  do i=1,nmax+1
 	  	do j=1,nmax+1
 	  	  do k=1,nmax+1
 	  	  	coeffs(i,j,k) = 0.d0
 	  	  	write(*,*) coeffs(i,j,k)
 	  	  end do
 	  	end do
 	  end do
 	  
 	  write(*,*) coeffs(1,3,28)
 	  coeffs(1,3,28) = 0.0d0
 	  write(*,*) coeffs(1,3,28)

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
 	      	    	if(i==1.AND.j==3.AND.k==28)then
 	      	    		if (kk==1) then
 	      	    		write(*,*) coeffs(i,j,k)
 	      	    		end if
 	      	    	end if
 	      	    end do
 	      	  end do
 	      	end do
 	      end do
 	    end do
 	  end do
 	  write(*,*) coeffs(1,3,28)
 	end subroutine
 	
  !----------------------------------------------------------------------
 	  ! RECONSTRUCT_CUBE reconstructs the data in the case of a nice cube.
 	  	  ! This follows strictly the algorithm and optimization layed out
 	  	  ! in Fluke et. al.
  !----------------------------------------------------------------------	
  	  subroutine reconstruct_cube(coeffs,Ng,x_max, recon)
  	  	integer, intent(in)	:: Ng
  	  	real(8), intent(in)	:: x_max
  	  	real(8), intent(in) :: coeffs(nmax+1,nmax+1,nmax+1)
  	  	
  	  	real(8), intent(out):: recon(Ng,Ng,Ng)
  	  	
  	  	integer	:: i,j,k,ii,jj,kk
  	  	real(8)	:: bases(Ng,nmax+1)
  	  	real(16)::x

  	  	recon(:,:,:) = 0.0d0
  	  	
  	  	do j = 1,Ng
  	  		x = -x_max+(j-0.5d0)*dx
  	  		call basis_vector(x,nmax,bases(j,:))
  	  	end do
  	  	
  	  	write(*,*) "IN RECONSTRUCT_CUBE"
  	  	write(*,*) 'Ng = ', Ng
  	  	write(*,*) 'x_max = ', x_max

  	  	do k=1,Ng
  	  		do j = 1,Ng
  	  			do i = 1,Ng
  	  				!write(*,*) "CELL:", i,j,k
  	  				do kk = 1,nmax+1
  	  					jj=1
  	  					ii=1
  	  					if (ii+jj+kk .GT. nmax+1)then
  	  						exit
  	  					end if
  	  					do jj = 1,nmax+1
  	  						ii=1
  	  						if (ii+jj+kk .GT. nmax+1)then
  	  							exit
  	  						end if	
  	  						do ii = 1,nmax+1
  	  						if (ii+jj+kk .GT. nmax+1)then
  	  							exit
  	  						end if
  	  							recon(i,j,k) = recon(i,j,k) + &
  	  								& coeffs(ii,jj,kk)*bases(i,ii)*&
  	  								&bases(j,jj)*bases(k,kk)
  	  							!if ((ii.LT.10).AND.jj==1.AND.kk==1)then	
  	  							!	write(*,*) recon(i,j,k)
  	  							!end if
  	  						end do
  	  					end do
  	  				end do
  	  			end do
  	  		end do
  	  	end do
  	  	
  	  	!write(*,*) 'reconstructed cube:'
  	    !write(*,*) recon
  	  end subroutine
  	  
end module
