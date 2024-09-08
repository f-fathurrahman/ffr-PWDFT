
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putepsinv(iq,epsi)
use modmain
implicit none
! arguments
integer, intent(in) :: iq
complex(8), intent(in) :: epsi(ngrf,ngrf,nwrf)
! local variables
integer recl,i
! determine the record length for EPSINV.OUT
inquire(iolength=recl) vql(:,iq),ngrf,nwrf,epsi
!$OMP CRITICAL(u180)
do i=1,2
  open(180,file='EPSINV.OUT',form='UNFORMATTED',access='DIRECT',recl=recl, &
   err=10)
  write(180,rec=iq,err=10) vql(:,iq),ngrf,nwrf,epsi
  close(180)
  exit
10 continue
  if (i.eq.2) then
    write(*,*)
    write(*,'("Error(putepsinv): unable to write to EPSINV.OUT")')
    write(*,*)
    stop
  end if
  close(180)
end do
!$OMP END CRITICAL(u180)
return
end subroutine

