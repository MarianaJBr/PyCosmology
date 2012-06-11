module f2pytest
	implicit none
	
	
	real, allocatable	:: a(:,:)
	
	contains
	
	subroutine allocate_a()
		
		allocate(a(5,5))
	end subroutine
	
	subroutine print_a
		write(*,*) a
	end subroutine
end module
	
		
