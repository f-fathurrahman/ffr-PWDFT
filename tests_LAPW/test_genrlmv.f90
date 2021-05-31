program test_genrlmv
  implicit none
  integer :: lmax, lmmax
  real(8) :: v(3)
  real(8), allocatable :: rlm(:)
  integer :: i

  lmax = 2
  lmmax = (lmax+1)**2
  allocate(rlm(lmmax))

  v = (/ 1.0d0, 1.1d0, 1.2d0 /)
  call genrlmv(lmax, v, rlm)

  do i=1,lmmax
    write(*,'(1x,I4,F18.10)') i, rlm(i)
  enddo

  deallocate(rlm)

end program