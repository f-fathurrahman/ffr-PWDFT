subroutine my_genpmat()
use modmain
implicit none
! local variables
integer :: ik

!ffr: I'm expecting that this will call my_genpmatk
do ik=1,nkpt
  write(*,'("Info(my_genpmat): ",I6," of ",I6," k-points")') ik,nkpt
  call my_putpmat(ik)
end do
return
end subroutine

