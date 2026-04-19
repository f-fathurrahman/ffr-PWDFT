subroutine my_phfext(iq,is,ia,ip,fext)
  use modmain, only: ivq, ngridq
  implicit none
  ! arguments
  integer, intent(in) :: iq,is,ia,ip
  character(*), intent(out) :: fext
  ! local variables
  integer i,j,m(3),n(3)
  ! external functions
  integer, external :: gcd
  !
  do i = 1,3
    if (ivq(i,iq) /= 0) then
      j = gcd(ivq(i,iq), ngridq(i))
      m(i) = ivq(i,iq)/j
      n(i) = ngridq(i)/j
    else
      m(i) = 0
      n(i) = 0
    endif
  enddo
  write(fext,'("_Q",2I2.2,"_",2I2.2,"_",2I2.2,"_S",I2.2,"_A",I3.3,"_P",I1,".OUT")') &
            & m(1),n(1), m(2),n(2), m(3),n(3), is, ia, ip
  return
end subroutine

