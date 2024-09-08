
! Copyright (C) 2019 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine getephmat(iq,ik,ephmat)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq,ik
complex(8), intent(out) :: ephmat(nstsv,nstsv,nbph)
! local variables
integer recl,n,i
integer nstsv_,nbph_
real(8) vql_(3),vkl_(3),t1
! find the record length
inquire(iolength=recl) vql_(:),vkl_(:),nstsv_,nbph_,ephmat
! record number
n=(iq-1)*nkptnr+ik
!$OMP CRITICAL(u240)
do i=1,2
  open(240,file='EPHMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl, &
   err=10)
  read(240,rec=n,err=10) vql_,vkl_,nstsv_,nbph_,ephmat
  exit
10 continue
  if (i.eq.2) then
    write(*,*)
    write(*,'("Error(getephmat): unable to read from EPHMAT.OUT")')
    write(*,*)
    stop
  end if
  close(240)
end do
!$OMP END CRITICAL(u240)
t1=abs(vql(1,iq)-vql_(1))+abs(vql(2,iq)-vql_(2))+abs(vql(3,iq)-vql_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getephmat): differing vectors for q-point ",I8)') iq
  write(*,'(" current    : ",3G18.10)') vql(:,iq)
  write(*,'(" EPHMAT.OUT : ",3G18.10)') vql_
  write(*,*)
  stop
end if
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getephmat): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EPHMAT.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getephmat): differing nstsv for q-point, k-point ",2I8)') &
   iq,ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" EPHMAT.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
if (nbph.ne.nbph_) then
  write(*,*)
  write(*,'("Error(getephmat): differing nbph for q-point, k-point ",2I8)') &
   iq,ik
  write(*,'(" current    : ",I8)') nbph
  write(*,'(" EPHMAT.OUT : ",I8)') nbph_
  write(*,*)
  stop
end if
return
end subroutine

