
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_sp_2b
! !INTERFACE:
subroutine ggamt_sp_2b(is,np,g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2, &
 dxdgd2,dxdgud,dcdgu2,dcdgd2,dcdgud)
! !USES:
use modmain
! !DESCRIPTION:
!   Post processing step of muffin-tin gradients for GGA type 2. See routine
!   {\tt ggamt\_sp\_2a} for full details.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is,np
real(8), intent(in) :: g2up(np),g2dn(np)
real(8), intent(in) :: gvup(np,3),gvdn(np,3)
real(8), intent(inout) :: vxup(np),vxdn(np)
real(8), intent(inout) :: vcup(np),vcdn(np)
real(8), intent(in) :: dxdgu2(np),dxdgd2(np),dxdgud(np)
real(8), intent(in) :: dcdgu2(np),dcdgd2(np),dcdgud(np)
! local variables
integer nr,nri,i
! allocatable arrays
real(8), allocatable :: rfmt1(:),rfmt2(:),grfmt(:,:)
allocate(rfmt1(np),rfmt2(np),grfmt(np,3))
nr=nrmt(is)
nri=nrmti(is)
!------------------!
!     exchange     !
!------------------!
! convert dxdgu2 to spherical harmonics
call rfsht(nr,nri,dxdgu2,rfmt1)
! compute grad dxdgu2
call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rfmt1,np,grfmt)
! (grad dxdgu2).(grad rhoup) in spherical coordinates
rfmt1(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt2)
  rfmt1(1:np)=rfmt1(1:np)+rfmt2(1:np)*gvup(1:np,i)
end do
vxup(1:np)=vxup(1:np)-2.d0*(rfmt1(1:np)+dxdgu2(1:np)*g2up(1:np)) &
 -dxdgud(1:np)*g2dn(1:np)
! convert dxdgd2 to spherical harmonics
call rfsht(nr,nri,dxdgd2,rfmt1)
! compute grad dxdgd2
call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rfmt1,np,grfmt)
! (grad dxdgd2).(grad rhodn) in spherical coordinates
rfmt1(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt2)
  rfmt1(1:np)=rfmt1(1:np)+rfmt2(1:np)*gvdn(1:np,i)
end do
vxdn(1:np)=vxdn(1:np)-2.d0*(rfmt1(1:np)+dxdgd2(1:np)*g2dn(1:np)) &
 -dxdgud(1:np)*g2up(1:np)
! convert dxdgud to spherical harmonics
call rfsht(nr,nri,dxdgud,rfmt1)
! compute grad dxdgud
call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rfmt1,np,grfmt)
! (grad dxdgud).(grad rhodn) and (grad dxdgud).(grad rhoup)
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt1)
  vxup(1:np)=vxup(1:np)-rfmt1(1:np)*gvdn(1:np,i)
  vxdn(1:np)=vxdn(1:np)-rfmt1(1:np)*gvup(1:np,i)
end do
!---------------------!
!     correlation     !
!---------------------!
! convert dcdgu2 to spherical harmonics
call rfsht(nr,nri,dcdgu2,rfmt1)
! compute grad dcdgu2
call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rfmt1,np,grfmt)
! (grad dcdgu2).(grad rhoup) in spherical coordinates
rfmt1(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt2)
  rfmt1(1:np)=rfmt1(1:np)+rfmt2(1:np)*gvup(1:np,i)
end do
vcup(1:np)=vcup(1:np)-2.d0*(rfmt1(1:np)+dcdgu2(1:np)*g2up(1:np)) &
 -dcdgud(1:np)*g2dn(1:np)
! convert dcdgd2 to spherical harmonics
call rfsht(nr,nri,dcdgd2,rfmt1)
! compute grad dcdgd2
call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rfmt1,np,grfmt)
! (grad dcdgd2).(grad rhodn) in spherical coordinates
rfmt1(1:np)=0.d0
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt2)
  rfmt1(1:np)=rfmt1(1:np)+rfmt2(1:np)*gvdn(1:np,i)
end do
vcdn(1:np)=vcdn(1:np)-2.d0*(rfmt1(1:np)+dcdgd2(1:np)*g2dn(1:np)) &
 -dcdgud(1:np)*g2up(1:np)
! convert dcdgud to spherical harmonics
call rfsht(nr,nri,dcdgud,rfmt1)
! compute grad dcdgud
call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rfmt1,np,grfmt)
! (grad dcdgud).(grad rhodn) and (grad dcdgud).(grad rhoup)
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),rfmt1)
  vcup(1:np)=vcup(1:np)-rfmt1(1:np)*gvdn(1:np,i)
  vcdn(1:np)=vcdn(1:np)-rfmt1(1:np)*gvup(1:np,i)
end do
deallocate(rfmt1,rfmt2,grfmt)
return
end subroutine
!EOC

