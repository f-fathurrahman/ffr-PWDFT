
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getoccsv(fext,ikp,vpl,occsvp)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ikp
real(8), intent(in) :: vpl(3)
real(8), intent(out) :: occsvp(nstsv)
! local variables
integer isym,ik
integer recl,nstsv_
real(8) vkl_(3),t1
if (ikp.gt.0) then
  ik=ikp
else
! find the k-point number
  call findkpt(vpl,isym,ik)
end if
! find the record length
inquire(iolength=recl) vkl_,nstsv_,occsvp
!$OMP CRITICAL(u130)
open(130,file='OCCSV'//trim(fext),form='UNFORMATTED',access='DIRECT', &
 action='READ',recl=recl)
read(130,rec=ik) vkl_,nstsv_,occsvp
close(130)
!$OMP END CRITICAL(u130)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getoccsv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" OCCSV.OUT  : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getoccsv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" OCCSV.OUT  : ",I8)') nstsv_
  write(*,*)
  stop
end if
return
end subroutine

