program test_sbessel
  implicit none
  integer, parameter :: lmax=4
  real(8) :: jl(0:lmax)
  real(8) :: x
  integer :: l
  
  x = 1.2d0
  CALL sbessel(lmax, x, jl)
  write(*,*) 'x = ', x
  do l=0,lmax
    write(*,'(1x,I4,F18.10)') l, jl(l)
  enddo
end program