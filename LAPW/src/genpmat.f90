
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genpmat
use modmain
implicit none
! local variables
integer :: ik

do ik=1,nkpt
  write(*,'("Info(genpmat): ",I6," of ",I6," k-points")') ik,nkpt
  call putpmat(ik)
end do
return
end subroutine

