module SimOperations                                             
  implicit none
  
  contains
  
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
    
    
    
      
end module

    
    
  
  
  
      
    
    
    
