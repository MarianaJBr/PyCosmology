
module shapelets
  implicit none
  
  !A module defining various subroutines and functions use in 3D shapelet
  !decomposition.
    
  	real(16), parameter	   ::pi= 3.141592653589  !Value of Pi
 	complex(16), parameter ::ci= (0.0d0,1.0d0)   !The Complex Unit (i)
 	real(16)				   :: beta_c
 	
 	contains

! 	!----------------------------------------------------------------------
! 	  ! HERMITE POLYNOMIAL
! 	!---------------------------------------------------------------------- 	
! 	subroutine Hermite(x,n,H)
! 	!Calculates the first n+1 polynomials (including 0). NOT USED IN THE
! 	!CUBE DECOMPOSITION
!
!    real(16), intent(in) ::x
!    integer, intent(in) ::n
!      
!    real(16),intent(out) ::H(n+1)
!    integer             ::i
!    
!    H(1) = 1
!    H(2) = 2.d0*x
!    
!    if (n.GT.2)then
!    do i=2,n
!      H(i+1) = 2.d0*x*H(i) - 2*i*H(i-1)
!    end do
!    end if
!
!  end subroutine
!    
!  !----------------------------------------------------------------------
! 	  ! PHI
! 	!----------------------------------------------------------------------	
! 	subroutine phi(x,n,p)
! 	!Calculates first n+1 phi's from Fluke et. al. NOT USED IN THE CUBE
! 	!DECOMPOSITION
! 	
!    real(16), intent(in) ::x
!    integer, intent(in) ::n
!    
!    real(16),intent(out) ::p(n+1)
!    real(16)             ::H(n+1)
!    integer             ::i
!    
!    call Hermite(x,n,H)
!    do i=0,n
!    
!    p(i+1) = (2**i * pi**0.5d0 * Factorial(i))**(-0.5d0) &
!                                        &* H(i+1)* exp(-x**2 /2)
!                      
!    end do
!    
!    end subroutine
  !----------------------------------------------------------------------
 	  ! BASIS_VECTOR - calculates bases iteratively
  !----------------------------------------------------------------------  
  	  subroutine basis_vector(x,n,bases)
  	  !Calculates the first n+1 basis functions at a point x, based on the
  	  !iterative method of Fluke et. al.
  	  	  
  	  	real(16), intent(in)	::x
  	  	integer, intent(in)	::n
  	  	
  	  	real(8), intent(out)	::bases(n+1)
  	  	integer				::i
  	  	
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
  	  
    
!  !----------------------------------------------------------------------
! 	  ! SHAPELET BASIS FUNCTIONS
!  !----------------------------------------------------------------------	  
! 	subroutine basis(x,n,b,Base)
! 	!Calculates the first n+1 basis functions at a point x, iterating on the
! 	!Hermite polynomials rather than the basis functions themselves. NOT USED
! 	!IN THE CUBE DECOMPOSITION.
! 		
! 	  real(16), intent(in) ::x
! 	  integer, intent(in) ::n
! 	  real(16), intent(in) ::b
! 	  real(16)			   :: Base(n+1),p(n+1)
! 	  
!    
! 	  call phi(x/b,n,p)
! 	  Base = b**(-1.d0/2.d0) * p
! 	 end subroutine

  !----------------------------------------------------------------------
 	  ! cube_ints - finds the integrals necessary for a cube
  !----------------------------------------------------------------------	  
 	subroutine cube_ints(x_max,Ng,n,Ints)
 	!Calculates the first n+1 Integral factors for all Ng cells. Beta must be
 	!specified beforehand. THIS IS NOT WORKING PROPERLY PROBABLY BECAUSE OF ERF.
 	
 	  integer, intent(in)  ::n,Ng
 	  real(16), intent(in)	::x_max
 	  
 	  integer              ::i,j
 	  real(8)              ::bases(Ng+1,n+1)
 	  real(8), intent(out) ::Ints(Ng,n+1)
 	  real(16)				:: x_a, x_b, erf_a, erf_b, a, b

 	  do i=1,Ng+1
 	  	a = -x_max + (i-1)*2.d0*x_max/Ng
 	  	call basis_vector(a,n,bases(i,:))
 	  end do
 	  
 	  do i=1,Ng
 	  	a = -x_max + (i-1)*2.d0*x_max/Ng
 	  	b = a + 2.d0*x_max/Ng 	 
 	  	x_a = a/(beta_c*sqrt(2.d0))
 	  	x_b = b/(beta_c*sqrt(2.d0))
 	  	
 	  	erf_a = erf(x_a)
 	  	erf_b = erf(x_b)
		write(*,'(A9,es28.20)') "erf(a) = ", log(erf_a)
		write(*,'(A9,es28.20)') "erf(b) = ", log(erf_b)
		
 	  	Ints(i,1) = (beta_c*pi**0.5d0 /2.d0)**(0.5d0) * (erf_b-erf_a)
		  
		Ints(i,2) = -sqrt(2.d0)*beta_c*(bases(i+1,1)-bases(i,1))
 	  
		do j=3,n+1
 	    	Ints(i,j) = -beta_c*sqrt(2.d0/n)*(bases(i+1,j-1)-bases(i,j-1))+&
 	                &sqrt((i-1)/real(i))*Ints(i,j-2)
 	    end do
 	  end do
 	  
 	  do i=1,n+1
 	  	write(*,'(es28.20)')Ints(1,i)
 	  end do
 	  	
 	end subroutine

  !----------------------------------------------------------------------
 	  ! COEFF_CUBE evaluates the coefficients in the case of a nice cube.
 	  	  ! This follows strictly the algorithm and optimization layed out
 	  	  ! in Fluke et. al.
  !----------------------------------------------------------------------	 	 
 	subroutine coeff_cube(f,Ng,x_max_8,cube)
 	  integer, intent(in)   :: Ng
 	  real(8), intent(in)   :: f(Ng,Ng,Ng)
 	  real(8), intent(in)   :: x_max_8
 	  
 	  
 	  real(8), intent(out)  :: cube((Ng-3)/2+1,(Ng-3)/2+1,(Ng-3)/2+1)
 	  
 	  real(8)				 ::Ints(Ng,(Ng-3)/2+1)
 	  real(8)               ::dx
 	  integer               ::i,j,k,nmax,ii,jj,kk
 	  real(16)				 :: x_max
 	  
 	  x_max = x_max_8
 	  dx = 2.d0*x_max/Ng
 	  nmax = (Ng-3)/2
 	  beta_c = x_max/sqrt(2.0d0*Ng)

 	  write(*,*) "x_max = ", x_max
 	  write(*,*) "Ng = ", Ng
 	  write(*,*) "nmax = ", nmax
 	  write(*,*) "dx = ", dx
 	  write(*,*) "beta = ", beta_c
 	  
 	  cube = 0.0d0
 	  call cube_ints(x_max,Ng,nmax,Ints)
 

 	  do k=1,nmax+1
 	  	i=1
 	  	j=1
 	  	!if (i+j+k .GT. nmax+1)then
 	    !  		exit
 	    !end if
 	    do j=1,nmax+1
 	    	i=1
 	    	!if (i+j+k .GT. nmax+1)then
 	       !		exit
 	      	!end if
 	      do i=1,nmax+1
 	      	!if (i+j+k .GT. nmax+1)then
 	      !		exit
 	      	!end if
 	      	do kk=1,Ng
 	      	  do jj = 1,Ng
 	      	    do ii = 1,Ng
 	      	    	cube(i,j,k) = cube(i,j,k) + &
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
  	  subroutine reconstruct_cube(cube,Ng,f,nmax,bounds)
  	  	integer, intent(in)	:: Ng,nmax
  	  	real(8), intent(in) :: f(nmax+1,nmax+1,nmax+1),bounds(2)
  	  	real(8), intent(out) :: cube(Ng,Ng,Ng)
  	  	
  	  	integer	:: i,j,k,ii,jj,kk
  	  	real(8)	:: bases(Ng,nmax+1)
  	  	real(16)::x
  	  	
  	  	
  	  	cube = 0.0d0
  	  	
  	  	do j = 1,Ng
  	  		x = bounds(1)+(j-0.5d0)*2.0*bounds(2)/Ng
  	  		call basis_vector(x,nmax,bases(j,:))
  	  	end do
  	  	
  	  	do k=1,Ng
  	  		do j = 1,Ng
  	  			do i = 1,Ng
  	  				do kk = 1,nmax+1
  	  					jj=1
  	  					ii=1
  	  					!if (ii+jj+kk .GT. nmax+1)then
  	  				!		exit
  	  				!	end if
  	  					do jj = 1,nmax+1
  	  						ii=1
  	  				!		if (ii+jj+kk .GT. nmax+1)then
  	  				!			exit
  	  				!		end if	
  	  						do ii = 1,nmax+1
  	  				!		if (ii+jj+kk .GT. nmax+1)then
  	  				!			exit
  	  				!		end if
  	  							cube(i,j,k) = cube(i,j,k) + &
  	  								& f(ii,jj,kk)*bases(i,ii)*&
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
  	  	real(16), intent(out)::cube_0(Ng,Ng,Ng), cube_1(Ng,Ng,Ng)
  	  	real(16), intent(out)::cube_2(Ng,Ng,Ng), cube_3(Ng,Ng,Ng)
  	  	
  	  	real(16)		:: x_max, x 
  	  	real(8)		:: base(5,Ng)
  	  	integer		::i,j,k
  	  	
  	  	
  	  	x_max = 3.d0
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
 	  
!	SUBROUTINE CERROR(Z,CER,acc)
!
!       ====================================================
!       Purpose: Compute error function erf(z) for a complex
!                argument (z=x+iy)
!       Input :  z   --- Complex argument
!       Output:  CER --- erf(z)
!       ====================================================
!
!	IMPLICIT COMPLEX(16)	:: (C,Z)
!	real(16) 				:: A0,PI
!	A0=CDABS(Z)
!	C0=CDEXP(-Z*Z)
!	PI=3.141592653589793D0
!	Z1=Z
!	IF (REAL(Z).LT.0.0) THEN
!	   Z1=-Z
!	ENDIF
!	IF (A0.LE.5.8D0) THEN    
!	   CS=Z1
!	   CR=Z1
!	   DO K=1,120
!	      CR=CR*Z1*Z1/(K+0.5D0)
!	      CS=CS+CR
!	      IF (CDABS(CR/CS).LT.1.0D-15)then
!	      	  GOTO 15
!	      end if
!	   end do
!15         CER=2.0D0*C0*CS/DSQRT(PI)
!	ELSE                              
!	   CL=1.0D0/Z1              
!	   CR=CL
!	   DO 20 K=1,13
!	      CR=-CR*(K-0.5D0)/(Z1*Z1)
!	      CL=CL+CR
!	      IF (CDABS(CR/CL).LT.1.0D-15)then
!	      	  GOTO 25
!	      end if
!	   end do
!25         CER=1.0D0-C0*CL/DSQRT(PI)
!	ENDIF
!	IF (REAL(Z).LT.0.0) THEN
!	   CER=-CER
!	ENDIF
!	RETURN
!	END
