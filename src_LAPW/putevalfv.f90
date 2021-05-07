
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevalfv(fext,ik,evalfv)
use modmain
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik
real(8), intent(in) :: evalfv(nstfv,nspnfv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstfv,nspnfv,evalfv
!$OMP CRITICAL(u120)
open(120,file='EVALFV'//trim(fext),form='UNFORMATTED',access='DIRECT', &
 action='WRITE',recl=recl)
write(120,rec=ik) vkl(:,ik),nstfv,nspnfv,evalfv
close(120)
!$OMP END CRITICAL(u120)
return
end subroutine

