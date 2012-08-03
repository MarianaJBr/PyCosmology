program WriteHealPix
	use pix_tools
	
	integer :: i,nside
	real(8)	:: vector(3)
	open(unit=1,file="healpix.dat")
	
	nside = 32
	write(1,'(i7)') 12*nside**2
	
	do i=0,12*nside**2-1
		call pix2vec_nest(nside, i, vector)
		write(1,'(3es20.8)') vector(1),vector(2),vector(3)
	end do
	close(1)
end program
