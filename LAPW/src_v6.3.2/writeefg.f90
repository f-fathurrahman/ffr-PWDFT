
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeefg
! !INTERFACE:
subroutine writeefg
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the electric field gradient (EFG) tensor for each atom, $\alpha$,
!   and writes it to the file {\tt EFG.OUT} along with its eigenvalues. The EFG
!   is defined by
!   $$ V^{\alpha}_{ij}\equiv\left.\frac{\partial^2 V'_{\rm C}({\bf r})}
!    {\partial{\bf r}_i\partial{\bf r}_j}\right|_{{\bf r}={\bf r}_{\alpha}}, $$
!   where $V'_{\rm C}$ is the Coulomb potential with the $l=m=0$ component
!   removed in each muffin-tin. The derivatives are computed explicitly using
!   the routine {\tt gradrfmt}.
!
! !REVISION HISTORY:
!   Created May 2004 (JKD)
!   Fixed serious problem, November 2006 (JKD)
!EOP
!BOC
implicit none
! local variables
integer, parameter :: lwork=10
integer is,ia,ias
integer nr,nri,ir
integer np,i,j,info
real(8) efg(3,3),a(3,3)
real(8) w(3),work(lwork)
! allocatable arrays
real(8), allocatable :: rfmt(:),grfmt1(:,:),grfmt2(:,:)
if (lmaxi.lt.2) then
  write(*,*)
  write(*,'("Error(writeefg): lmaxi too small for calculating the EFG : ",&
   &I4)') lmaxi
  write(*,'(" Run the ground-state calculation again with lmaxi >= 2")')
  write(*,*)
  stop
end if
! initialise universal variables
call init0
! read density and potentials from file
call readstate
! allocate local arrays
allocate(rfmt(npmtmax),grfmt1(npmtmax,3),grfmt2(npmtmax,3))
open(50,file='EFG.OUT',form='FORMATTED')
write(50,*)
write(50,'("(electric field gradient tensor is in Cartesian coordinates)")')
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
! remove the l=m=0 part of the potential
    rfmt(1:np)=vclmt(1:np,ias)
    i=1
    do ir=1,nri
      rfmt(i)=0.d0
      i=i+lmmaxi
    end do
    do ir=nri+1,nr
      rfmt(i)=0.d0
      i=i+lmmaxo
    end do
! compute the gradient of the Coulomb potential
    call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rfmt,npmtmax,grfmt1)
    do i=1,3
! compute the gradient of the gradient
      call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),grfmt1(:,i),npmtmax, &
       grfmt2)
      do j=1,3
        efg(i,j)=grfmt2(1,j)*y00
      end do
    end do
! symmetrise the EFG
    do i=1,3
      do j=i+1,3
        efg(i,j)=0.5d0*(efg(i,j)+efg(j,i))
        efg(j,i)=efg(i,j)
      end do
    end do
    write(50,*)
    write(50,'(" EFG tensor :")')
    do i=1,3
      write(50,'(3G18.10)') (efg(i,j),j=1,3)
    end do
    write(50,'(" trace : ",G18.10)') efg(1,1)+efg(2,2)+efg(3,3)
! diagonalise the EFG
    a(:,:)=efg(:,:)
    call dsyev('N','U',3,a,3,w,work,lwork,info)
    write(50,'(" eigenvalues :")')
    write(50,'(3G18.10)') w
  end do
end do
close(50)
write(*,*)
write(*,'("Info(writeefg): electric field gradient written to EFG.OUT")')
deallocate(rfmt,grfmt1,grfmt2)
return
end subroutine
!EOC

