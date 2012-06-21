
module shapelets
  implicit none
  
  !A module defining various subroutines and functions use in 3D shapelet
  !decomposition.
    
  	real(8), parameter	   ::pi= 3.141592653589  !Value of Pi
 	complex(8), parameter ::ci= (0.0d0,1.0d0)   !The Complex Unit (i)
 	real(8)				   :: beta
 	
 	contains
 	
 	!----------------------------------------------------------------------
 	  ! FACTORIAL
 	!----------------------------------------------------------------------
 	function factorial(n)
 	  integer, intent(in) ::n
 	  
 	  integer(8)          ::factorial
 	  integer             ::i
 	  
 	  factorial = 1
 	  if (n.EQ.0 .OR. n.EQ.1) then
 	     return
 	  else
 	
 	  
          do i=1,n
            factorial = factorial*i
          end do
      end if
 	end function

 	!----------------------------------------------------------------------
 	  ! HERMITE POLYNOMIAL
 	!---------------------------------------------------------------------- 	
 	subroutine Hermite(x,n,H)

    real(8), intent(in) ::x
    integer, intent(in) ::n
      
    real(8),intent(out) ::H(n+1)
    integer             ::i
    
    H(1) = 1
    H(2) = 2.d0*x
    
    if (n.GT.2)then
    do i=2,n
      H(i+1) = 2.d0*x*H(i) - 2*i*H(i-1)
    end do
    end if

  end subroutine
    
  !----------------------------------------------------------------------
 	  ! PHI
 	!----------------------------------------------------------------------	
 	subroutine phi(x,n,p)
 	
    real(8), intent(in) ::x
    integer, intent(in) ::n
    
    real(8),intent(out) ::p(n+1)
    real(8)             ::H(n+1)
    integer             ::i
    
    call Hermite(x,n,H)
    do i=0,n
    
    p(i+1) = (2**i * pi**0.5d0 * Factorial(i))**(-0.5d0) &
                                        &* H(i+1)* exp(-x**2 /2)
                      
    end do
    
    end subroutine
  !----------------------------------------------------------------------
 	  ! BASIS_VECTOR - calculates bases iteratively
  !----------------------------------------------------------------------  
  	  subroutine basis_vector(x,n,bases)
  	  	real(8), intent(in)	::x
  	  	integer, intent(in)	::n
  	  	
  	  	real(8), intent(out)	::bases(n+1)
  	  	integer				::i
  	  	
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
 	  ! SHAPELET BASIS FUNCTIONS
  !----------------------------------------------------------------------	  
 	subroutine basis(x,n,b,Base)
 	
 	  real(8), intent(in) ::x
 	  integer, intent(in) ::n
 	  real(8), intent(in) ::b
 	  real(8)			   :: Base(n+1),p(n+1)
 	  
    
 	  call phi(x/b,n,p)
 	  Base = b**(-1.d0/2.d0) * p
 	end subroutine
  
  !----------------------------------------------------------------------
 	  ! Integrals
 	!----------------------------------------------------------------------	  
 	subroutine Integral(a,b,beta,n,Ints)
 	
 	  real(8), intent(in)  ::a,b,beta
 	  integer, intent(in)  ::n
 	  
 	  integer              ::i
 	  real(8)              ::a_basis(n+1),b_basis(n+1)
 	  real(8), intent(out) ::Ints(n+1)
 	  
 	  call basis(a,n,beta,a_basis)
 	  call basis(b,n,beta,b_basis)
 	  
 	  Ints(1) = (beta*pi**0.5d0 /2.d0)**(0.5d0) * (erf(b/(beta*sqrt(2.d0)))&
 	  &-erf(a/(beta*sqrt(2.d0))))
 	  
 	  Ints(2) = -sqrt(2.d0)*beta*(b_basis(1)-a_basis(1))
 	  
 	  do i=3,n+1
 	    Ints(i) = -beta*sqrt(2.d0/n)*(b_basis(i-1)-a_basis(i-1))+&
 	                &sqrt((i-1)/real(i))*Ints(i-2)
 	  end do
 	end subroutine

  !----------------------------------------------------------------------
 	  ! cube_ints - finds the integrals necessary for a cube
 	!----------------------------------------------------------------------	  
 	subroutine cube_ints(x_max,Ng,n,Ints)
 	
 	  integer, intent(in)  ::n,Ng
 	  real(8), intent(in)	::x_max
 	  
 	  integer              ::i,j
 	  real(8)              ::bases(Ng+1,n+1),a,b
 	  real(8), intent(out) ::Ints(Ng,n+1)
 	  
 	  do i=1,Ng+1
 	  	a = -x_max + (i-1)*2.d0*x_max/Ng
 	  	call basis_vector(a,n,bases(i,:))
 	  end do
 	  
 	  do i=1,Ng
 	  	a = -x_max + (i-1)*2.d0*x_max/Ng
 	  	b = a + 2.d0*x_max/Ng 	  
		  
 	  	Ints(i,1) = (beta*pi**0.5d0 /2.d0)**(0.5d0) * (erf(b/(beta*sqrt(2.d0)))&
		&-erf(a/(beta*sqrt(2.d0))))
		  
		Ints(i,2) = -sqrt(2.d0)*beta*(bases(i+1,1)-bases(i,1))
 	  
		do j=3,n+1
 	    	Ints(i,j) = -beta*sqrt(2.d0/n)*(bases(i+1,j-1)-bases(i+1,j-1))+&
 	                &sqrt((i-1)/real(i))*Ints(i,j-2)
 	    end do
 	  end do
 	end subroutine
  
  !----------------------------------------------------------------------
 	  ! COEFF
  !----------------------------------------------------------------------	 	 
 	!function coeff(n,Ng,beta,f,bounds)
 	!
 	!  real(8), intent(in)   :: f(:,:,:)
 	!  real(8), intent(in)   :: bounds(6), beta(3)
 	!  integer, intent(in)   :: n(3), Ng(3)
 	!  
 	!  real(8)               ::coeff
 	!  
 	!  real(8)               ::dx,dy,dz, ax, ay, az, bx,by,bz
 	!  
 	!  integer               ::i,j,k
 	!  
 	!  dx = (bounds(2)-bounds(1))/Ng(1)
 	!  dy = (bounds(4)-bounds(3))/Ng(2)
 	!  dz = (bounds(6)-bounds(5))/Ng(3)
 	!  
 	!  do k=1,Ng(3)
 	!    az = bounds(5) + (i-1)*dz
 	!    bz = az + dz
 	!    do j=1,Ng(2)
 	!      ay = bounds(3) + (i-1)*dy
 	!      by = ay + dy
 	!      do i=1,Ng(1)
 	!        ax = bounds(1) + (i-1)*dx
 	!        bx = ax + dx
 	!        coeff = coeff + f(i,j,k)*Integral(ax,bx,beta(1),n(1))*&
 	!                &Integral(ay,by,beta(2),n(2))*Integral(az,bz,beta(3),n(3))
 	!      end do
 	!    end do
 	!  end do
 	!end function

  !----------------------------------------------------------------------
 	  ! COEFF_CUBE evaluates the coefficients in the case of a nice cube.
 	  	  ! This follows strictly the algorithm and optimization layed out
 	  	  ! in Fluke et. al.
  !----------------------------------------------------------------------	 	 
 	subroutine coeff_cube(f,Ng,x_max,cube)
 	  integer, intent(in)   :: Ng
 	  real(8), intent(in)   :: f(Ng,Ng,Ng)
 	  real(8), intent(in)   :: x_max
 	  
 	  
 	  real(8), intent(out)  :: cube((Ng-3)/2+1,(Ng-3)/2+1,(Ng-3)/2+1)
 	  
 	  real(8)				 ::Ints(Ng,(Ng-3)/2+1)
 	  real(8)               ::dx
 	  integer               ::i,j,k,nmax,ii,jj,kk


 	  dx = 2.d0*x_max/Ng
 	  nmax = (Ng-3)/2
 	  beta = x_max/sqrt(2.0d0*Ng)

 	  write(*,*) "x_max = ", x_max
 	  write(*,*) "Ng = ", Ng
 	  write(*,*) "nmax = ", nmax
 	  write(*,*) "dx = ", dx
 	  write(*,*) "beta = ", beta
 	  
 	  cube = 0.0d0
 	  call cube_ints(x_max,Ng,nmax,Ints)
 

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
  	  	real(8)	:: bases(Ng,nmax+1),x
  	  	
  	  	
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
 	  ! TEST BASIS_VECTOR - calculates bases iteratively
  !----------------------------------------------------------------------  
  	  subroutine test(Ng,cube_0,cube_1,cube_2,cube_3)
  	  	!This subroutine produces four cubes, based on the basis vectors
  	  	!with beta = 1, n=(0,0,0), (0,1,0),(2,0,2),(1,2,4).
  	  		
  	  	integer, intent(in)	::Ng
  	  	real(8), intent(out)::cube_0(Ng,Ng,Ng), cube_1(Ng,Ng,Ng)
  	  	real(8), intent(out)::cube_2(Ng,Ng,Ng), cube_3(Ng,Ng,Ng)
  	  	
  	  	real(8)		:: x_max, x 
  	  	real(8)		:: base(5,Ng)
  	  	integer		::i,j,k
  	  	
  	  	
  	  	x_max = 3.d0
  	  	beta = 1.d0
  	  	
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
 	  
 	  
 	  
 	  
 	  
 	  
 	  
 	  
 	  
 	  
 	  
 	  
