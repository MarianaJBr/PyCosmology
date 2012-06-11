include "SimOperations.f95"

program matrices
  use SimOperations
  implicit none
  
  real  ::  a(3,3), b(3,3), c(3,3),d(3)
  integer:: i,j
  
  a = reshape((/1,2,3,4,5,6,7,8,9/),(/3,3/))
  b = reshape((/1,2,3,4,5,6,7,8,9/),(/3,3/))
  
  do i=1,3
    c(:,i) = a(:,i) -b(i,:)
  end do
  
  do i=1,3
     write(*,*)c(i,:)
  end do
  
  d = (/2,9,5/)
  b = reshape((/1,1,1,2,2,2,3,3,3/),(/3,3/))
  
  call EigOrder(d,b)
  
  do i=1,3
    write(*,*) d(i)
  end do
  do i=1,3
    write(*,*) b(i,:)
  end do
  
  a = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))

  b = Inertia(1.0, a)
  
  do i=1,3
    write(*,*) b(i,:)
  end do
end program
