
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: olpistl
! !INTERFACE:
subroutine olpistl(ngp,igpig,ld,o)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   ld    : leading dimension of o (in,integer)
!   o     : overlap matrix (inout,complex(*))
! !DESCRIPTION:
!   Computes the interstitial contribution to the overlap matrix for the APW
!   basis functions. The overlap is given by
!   $$ O^{\rm I}({\bf G+k,G'+k})=\tilde{\Theta}({\bf G-G'}), $$
!   where $\tilde{\Theta}$ is the characteristic function. See routine
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp,igpig(ngkmax)
integer, intent(in) :: ld
complex(8), intent(inout) :: o(*)
! local variables
integer i1,i2,i3,j1,j2,j3
integer ig,i,j,k
do j=1,ngp
  k=(j-1)*ld
  ig=igpig(j)
  j1=ivg(1,ig); j2=ivg(2,ig); j3=ivg(3,ig)
  do i=1,j
    k=k+1
    ig=igpig(i)
    i1=ivg(1,ig)-j1; i2=ivg(2,ig)-j2; i3=ivg(3,ig)-j3
    o(k)=o(k)+cfunig(ivgig(i1,i2,i3))
  end do
end do
return
end subroutine
!EOC

