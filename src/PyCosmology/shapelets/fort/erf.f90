FUNCTION my_erf(x)
	REAL(8)				:: my_erf
	real(8), intent(in)	:: x	
	REAL(8)				:: gammp
	
	if(x.lt.0.)then
	  my_erf=-gammp(.5d0,x**2)
	else
	  my_erf=gammp(.5d0,x**2)
	endif

end function


FUNCTION gammp(a,x)
	REAL(8), intent(in)	:: a,x
	real(8)				:: gammp
	REAL(8)				:: gammcf,gamser,gln
	
	if(x.lt.0..or.a.le.0.)then
		write(*,*) 'bad arguments in gammp'
		stop
	end if
	if(x.lt.a+1.)then
	  call gser(gamser,a,x,gln)
	  gammp=gamser
	else
	  call gcf(gammcf,a,x,gln)
	  gammp=1.-gammcf
	endif

end function



SUBROUTINE gcf(gammcf,a,x,gln)
	INTEGER, parameter 	:: ITMAX = 500
	REAL(8), intent(in)	:: a,x
	real(8), parameter	:: EPS = 3.0E-22,FPMIN=1.0E-30
	real(8), intent(out)	:: gammcf,gln
	
	INTEGER 	:: i
	REAL(8)	:: an,b,c,d,del,h,gammln
	gln=gammln(a)
	b=x+1.-a
	c=1./FPMIN
	d=1./b
	h=d
	do i=1,ITMAX
	  an=-i*(i-a)
	  b=b+2.
	  d=an*d+b
	  if(abs(d).lt.FPMIN)d=FPMIN
	  c=b+an/c
	  if(abs(c).lt.FPMIN)c=FPMIN
	  d=1./d
	  del=d*c
	  h=h*del
	  if(abs(del-1.).lt.EPS)goto 1
	end do
	
	write(*,*) "error in gcf (a too large or ITMAX too small)"
	stop
		  
	1     gammcf=exp(-x+a*log(x)-gln)*h

end subroutine
      
SUBROUTINE gser(gamser,a,x,gln)
	INTEGER, parameter 	:: ITMAX = 500
	REAL(8), intent(in) 	:: a,x
	real(8), intent(out)	:: gamser,gln
	real(8), parameter	:: EPS = 3.0E-22
	
	INTEGER 	:: n
	REAL(8)	:: ap,del,sum,gammln
	
	gln=gammln(a)
	if(x.le.0.)then
	  if(x.lt.0.)then
		write(*,*) 'x < 0 in gser'
		stop
	  end if
	  gamser=0.
	  return
	endif
	ap=a
	sum=1./a
	del=sum
	do n=1,ITMAX
	  ap=ap+1.
	  del=del*x/ap
	  sum=sum+del
	  if(abs(del).lt.abs(sum)*EPS)goto 1
	end do
	write(*,*) 'a too large, ITMAX too small in gser'
	stop
	1     gamser=sum*exp(-x+a*log(x)-gln)
end subroutine

      
FUNCTION gammln(xx)
	REAL(8), intent(in)	:: xx
	real(8)				:: gammln
	INTEGER 			:: j
	REAL(8) 			:: ser,tmp,x,y
	real(8), save		:: stp, cof(6)
	stp = 2.5066282746310005d0
	cof = (/76.18009172947146d0,-86.50532032941677d0,&
	&24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
	&-.5395239384953d-5/)
	
	x=xx
	y=x
	tmp=x+5.5d0
	tmp=(x+0.5d0)*log(tmp)-tmp
	ser=1.000000000190015d0
	do j=1,6
	  y=y+1.d0
	  ser=ser+cof(j)/y
	end do
	gammln=tmp+log(stp*ser/x)
end function