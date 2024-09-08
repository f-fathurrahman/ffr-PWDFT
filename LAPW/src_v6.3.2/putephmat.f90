
! Copyright (C) 2019 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine putephmat(iq,ik,ephmat)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq,ik
complex(8), intent(in) :: ephmat(nstsv,nstsv,nbph)
! local variables
integer recl,n,i
! determine the record length
inquire(iolength=recl) vql(:,1),vkl(:,1),nstsv,nbph,ephmat
! record number
n=(iq-1)*nkptnr+ik
!$OMP CRITICAL(u240)
do i=1,2
  open(240,file='EPHMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl, &
   err=10)
  write(240,rec=n,err=10) vql(:,iq),vkl(:,ik),nstsv,nbph,ephmat
  close(240)
  exit
10 continue
  if (i.eq.2) then
    write(*,*)
    write(*,'("Error(putephmat): unable to write to EPHMAT.OUT")')
    write(*,*)
    stop
  end if
  close(240)
end do
!$OMP END CRITICAL(u240)
return
end subroutine

