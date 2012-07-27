module AuxiliaryOperations
 implicit none
 
 contains

 
    
    
    

!----------------------------------------------------------------------
 	  ! ReCentre - Recentres around the main particle of a group.
  !---------------------------------------------------------------------- 
  subroutine ReCentre(pos, BoxSize, centre)
    real(8), intent(in)    :: BoxSize
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

  !----------------------------------------------------------------------
 	  ! Inertia - Calculates the Moment Of Inertia Tensor.
  !----------------------------------------------------------------------  
  
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

  !----------------------------------------------------------------------
 	  ! EigOrder - Orders a set of Eigenvalues and vectors according to 
 	  	  !			the size of the Eienvalues.
  !----------------------------------------------------------------------  
  
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

  !----------------------------------------------------------------------
 	  ! GenRotate - Rotates 3D data in a general rotation.
  !----------------------------------------------------------------------  
  subroutine GenRotate(pos,u,theta,pos_rot)
    real, intent(in)    ::u(:), theta,pos(:,:)
    real, intent(out) :: pos_rot(:,:)
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
      pos_rot(:,i) = matmul(ROT,pos(:,i))
    end do
  end subroutine

  !----------------------------------------------------------------------
 	  ! LogBins - Bins data in log bins.
  !----------------------------------------------------------------------
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
    
    

  !----------------------------------------------------------------------
 	  ! jacobi - finds the Jacobian of a matrix.
  !----------------------------------------------------------------------
    
    
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
  
  !----------------------------------------------------------------------
 	  ! sort - sorts arr in place in ascending order.
  !----------------------------------------------------------------------
      
      SUBROUTINE sort2(n,arr,brr)
      INTEGER n,M,NSTACK
      REAL arr(n),brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,b,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          temp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
end module
