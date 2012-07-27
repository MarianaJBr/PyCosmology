include 'shapelets.f90'

program shapelet_test
 use shapelets
 implicit none
 
 real(8)	:: reconstruction_int(51,51,51)
 real(8)	::coeffs(25,25,25)
 real(8)	:: recon(51,51,51)
 
 call cyclic_test((/0,0,0/),1,1.3d0,51,reconstruction_int,coeffs,recon)
 
end program
