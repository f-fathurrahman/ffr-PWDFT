subroutine info_lattice()
  USE m_lattice
  USE m_constants, ONLY: pi
  implicit none
  integer :: i
  real(8) :: m(3,3)

  write(*,*) 'avec matrix'
  do i=1,3
    write(*,'(1x,3F18.10)') avec(i,:)
  enddo

  m = matmul(avec, transpose(bvec))/(2.d0*pi)
  write(*,*) 'avec*bvec matrix = '
  do i=1,3
    write(*,'(1x,3F18.10)') m(i,:)
  enddo


end subroutine