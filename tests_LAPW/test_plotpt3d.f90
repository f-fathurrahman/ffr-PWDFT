program main
  IMPLICIT NONE
  real(8) :: vclp3d(3,0:3)
  integer :: np3d(3)
  real(8), allocatable :: vpl(:,:)

  vclp3d(:,:) = 0.d0
  vclp3d(1,1) = 1.d0
  vclp3d(2,2) = 1.d0
  vclp3d(3,3) = 1.d0

  np3d(:) = (/ 3,3,3 /)
  
  allocate( vpl(3,np3d(1)*np3d(2)*np3d(3)) )

  call my_plotpt3d(vclp3d, np3d, vpl)

  deallocate( vpl )
end program

include 'my_plotpt3d.f90'