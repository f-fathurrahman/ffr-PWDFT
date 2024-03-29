include 'rf_interp.f90'

subroutine my_func(x, res)
  implicit none 
  real(8), intent(in) :: x
  real(8), intent(out) :: res
  res = exp(-0.1*x)
  return
end subroutine

program test_rf_interp
  implicit none 
  ! arguments
  INTEGER, parameter :: ni=10
  REAL(8) :: xi(ni), fi(ni)
  INTEGER, parameter :: no=7
  REAL(8) :: xo(no), fo(no)
  integer :: i
  real(8) :: dx, L, fxo

  L = 1.d0
  dx = L/dble(ni-1)
  do i=1,ni
    xi(i) = (i-1)*dx
    call my_func(xi(i), fi(i))
    write(*,'(1x,I4,2F18.10)') i, xi(i), fi(i)
  enddo

  dx = L/dble(no-1)
  do i=1,no
    xo(i) = (i-1)*dx
    !write(*,'(1x,I4,F18.10)') i, xo(i)
  enddo

  call rf_interp(ni, xi, fi, no, xo, fo)

  do i=1,no
    call my_func(xo(i), fxo)
    write(*,'(1x,I4,3F18.10)') i, xo(i), fo(i), fxo
  enddo

end program