
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writegclq
! !INTERFACE:
subroutine writegclq
! !USES:
use modmain
! !DESCRIPTION:
!   Outputs the volume-averaged integral of $4\pi/q^2$ in the small
!   parallelepiped around each discrete $q$-point to the file {\tt GCLQ.OUT}.
!   These represent the regularised Coulomb Green's function in reciprocal
!   space for small $q$. See the routine gengclq.
!
! !REVISION HISTORY:
!   Created June 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer iq
open(50,file='GCLQ'//trim(filext),form='FORMATTED')
write(50,'(I6," : nqpt; q-point, vql, gclq below")') nqpt
do iq=1,nqpt
  write(50,'(I6,4G18.10)') iq,vql(:,iq),gclq(iq)
end do
close(50)
return
end subroutine
!EOC

