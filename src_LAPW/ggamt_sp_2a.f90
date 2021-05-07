
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_sp_2a
! !INTERFACE:
subroutine ggamt_sp_2a(is,np,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the muffin-tin gradients $\nabla^2\rho^{\uparrow}$,
!   $\nabla^2\rho^{\downarrow}$, $\nabla\rho^{\uparrow}$,
!   $\nabla\rho^{\downarrow}$, $(\nabla\rho^{\uparrow})^2$,
!   $(\nabla\rho^{\downarrow})^2$ and
!   $\nabla\rho^{\uparrow}\cdot\nabla\rho^{\downarrow}$, which are passed in to
!   GGA functional subroutines of type 2. The exchange-correlation energy in
!   these routines has the functional form
!   $$ E_{xc}[\rho^{\uparrow},\rho^{\downarrow}]=\int d^3r\,\hat{\epsilon}_{xc}
!    \bigl(\rho^{\uparrow}({\bf r}),\rho^{\downarrow}({\bf r}),
!    (\nabla\rho^{\uparrow}({\bf r}))^2,(\nabla\rho^{\downarrow}({\bf r}))^2,
!    \nabla\rho^{\uparrow}({\bf r})
!    \cdot\nabla\rho^{\downarrow}({\bf r})\bigr), $$
!   where $\hat{\epsilon}_{xc}({\bf r})=\epsilon_{xc}({\bf r})\rho({\bf r})$ is
!   the xc energy per unit volume, with $\epsilon_{xc}$ being the xc energy per
!   electron, and $\rho=\rho^{\uparrow}+\rho^{\downarrow}$. From the gradients
!   above, type 2 GGA routines return $\epsilon_{xc}$, but not directly the xc
!   potentials. Instead they generate the derivatives
!   $\partial\hat{\epsilon}_{xc}/\partial\rho^{\uparrow}({\bf r})$,
!   $\partial\hat{\epsilon}_{xc}/\partial(\nabla\rho^{\uparrow}({\bf r}))^2$,
!   and the same for down spin, as well as
!   $\partial\hat{\epsilon}_{xc}/\partial(\nabla\rho^{\uparrow}({\bf r})
!   \cdot\nabla\rho^{\downarrow}({\bf r}))$. In a post-processing step invoked
!   by {\tt ggamt\_sp\_2b}, integration by parts is used to obtain the xc
!   potential explicitly with
!   \begin{align*}
!    V_{xc}^{\uparrow}({\bf r})=&\frac{\partial\hat{\epsilon}_{xc}}{\partial
!    \rho^{\uparrow}({\bf r})}-2\left(\nabla\frac{\partial\hat{\epsilon}_{xc}}
!    {\partial(\nabla\rho^{\uparrow})^2}\right)\cdot\nabla\rho^{\uparrow}
!    -2\frac{\hat{\epsilon}_{xc}}{\partial(\nabla\rho^{\uparrow})^2}\nabla^2
!    \rho^{\uparrow}\\
!    &-\left(\nabla\frac{\hat{\epsilon}_{xc}}{\partial(\nabla\rho^{\uparrow}
!    \cdot\nabla\rho^{\downarrow})}\right)\cdot\nabla\rho^{\downarrow}
!    -\frac{\partial\hat{\epsilon}_{xc}}{\partial(\nabla\rho^{\uparrow}\cdot
!    \nabla\rho^{\downarrow})}\nabla^2\rho^{\downarrow},
!   \end{align*}
!   and similarly for $V_{xc}^{\downarrow}$.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is,np
real(8), intent(in) :: rhoup(np),rhodn(np)
real(8), intent(out) :: g2up(np),g2dn(np)
real(8), intent(out) :: gvup(np,3),gvdn(np,3)
real(8), intent(out) :: gup2(np),gdn2(np),gupdn(np)
! local variables
integer nr,nri,i
! allocatable arrays
real(8), allocatable :: rfmt1(:),rfmt2(:),grfmt(:,:)
allocate(rfmt1(np),rfmt2(np),grfmt(np,3))
nr=nrmt(is)
nri=nrmti(is)
!----------------!
!     rho up     !
!----------------!
! convert rhoup to spherical harmonics
call rfsht(nr,nri,rhoup,rfmt1)
! compute grad^2 rhoup in spherical coordinates
call grad2rfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rlmt(:,-2,is),rfmt1,rfmt2)
call rbsht(nr,nri,rfmt2,g2up)
! grad rhoup in spherical coordinates
call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rfmt1,np,grfmt)
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),gvup(:,i))
end do
! (grad rhoup)^2
gup2(1:np)=gvup(1:np,1)**2+gvup(1:np,2)**2+gvup(1:np,3)**2
!------------------!
!     rho down     !
!------------------!
! convert rhodn to spherical harmonics
call rfsht(nr,nri,rhodn,rfmt1)
! compute grad^2 rhodn in spherical coordinates
call grad2rfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rlmt(:,-2,is),rfmt1,rfmt2)
call rbsht(nr,nri,rfmt2,g2dn)
! grad rhodn in spherical coordinates
call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rfmt1,np,grfmt)
do i=1,3
  call rbsht(nr,nri,grfmt(:,i),gvdn(:,i))
end do
! (grad rhodn)^2
gdn2(1:np)=gvdn(1:np,1)**2+gvdn(1:np,2)**2+gvdn(1:np,3)**2
! (grad rhoup).(grad rhodn)
gupdn(1:np)=gvup(1:np,1)*gvdn(1:np,1) &
           +gvup(1:np,2)*gvdn(1:np,2) &
           +gvup(1:np,3)*gvdn(1:np,3)
deallocate(rfmt1,rfmt2,grfmt)
return
end subroutine
!EOC

