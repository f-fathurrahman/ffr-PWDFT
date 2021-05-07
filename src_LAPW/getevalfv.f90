
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevalfv(fext,ikp,vpl,evalfv)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ikp
real(8), intent(in) :: vpl(3)
real(8), intent(out) :: evalfv(nstfv,nspnfv)
! local variables
integer isym,ik
integer recl,nstfv_,nspnfv_
real(8) vkl_(3),t1
if (ikp.gt.0) then
  ik=ikp
else
! find the k-point number
  call findkpt(vpl,isym,ik)
end if
! find the record length
inquire(iolength=recl) vkl_,nstfv_,nspnfv_,evalfv
!$OMP CRITICAL(u120)
open(120,file='EVALFV'//trim(fext),form='UNFORMATTED',access='DIRECT', &
 action='READ',recl=recl)
read(120,rec=ik) vkl_,nstfv_,nspnfv_,evalfv
close(120)
!$OMP END CRITICAL(u120)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevalfv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVALFV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstfv.ne.nstfv_) then
  write(*,*)
  write(*,'("Error(getevalfv): differing nstfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstfv
  write(*,'(" EVALFV.OUT : ",I8)') nstfv_
  write(*,*)
  stop
end if
if (nspnfv.ne.nspnfv_) then
  write(*,*)
  write(*,'("Error(getevalfv): differing nspnfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nspnfv
  write(*,'(" EVALFV.OUT : ",I8)') nspnfv_
  write(*,*)
  stop
end if
return
end subroutine

