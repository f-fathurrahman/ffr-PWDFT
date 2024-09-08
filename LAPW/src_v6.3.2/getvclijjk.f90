
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getvclijjk
! !INTERFACE:
subroutine getvclijjk(ikp,vclijjk)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced set (in,integer)
!   vclijjk : Coulomb matrix elements (out,complex(nstsv,nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Retrieves Coulomb matrix elements of the type $(i-jj-k)$ from the file
!   {\tt VCLIJJK.OUT}.
!
! !REVISION HISTORY:
!   Created 2009 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(out) :: vclijjk(nstsv,nstsv,nstsv,nkpt)
! local variables
integer recl,ios
! determine record length
inquire(iolength=recl) vclijjk
!$OMP CRITICAL(u262)
open(262,file='VCLIJJK.OUT',form='UNFORMATTED',access='DIRECT',recl=recl, &
 iostat=ios)
if (ios.ne.0) then
  write(*,*)
  write(*,'("Error(getvclijjk): error opening file VCLIJJK.OUT")')
  write(*,*)
  stop
end if
read(262,rec=ikp) vclijjk
close(262)
!$OMP END CRITICAL(u262)
return
end subroutine
!EOC

