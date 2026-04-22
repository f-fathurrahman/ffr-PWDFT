
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlxbse
use modmain
implicit none
! local variables
integer ik2
do ik2=1,nkptnr
  write(*,'("Info(hmlxbse): ",I6," of ",I6," k-points")') ik2,nkptnr
  call hmlxbsek(ik2)
end do
return
end subroutine

