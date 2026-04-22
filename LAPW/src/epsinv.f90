
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine epsinv
use modmain
implicit none
! local variables
integer iq,ik,ig,iw
integer n,nthd
! allocatable arrays
integer(8), allocatable :: lock(:)
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:),epsi(:,:,:)

! allocate local arrays
allocate(vgqc(3,ngrf),gqc(ngrf),gclgq(ngrf))
allocate(jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
allocate(epsi(ngrf,ngrf,nwrf))

write(*,*)

! loop over q-points
do iq=1,nqpt
  write(*,'("Info(epsinv): ",I6," of ",I6," q-points")') iq,nqpt
! generate the G+q-vectors and related quantities
  call gengqrf(vqc(:,iq),vgqc,gqc,jlgqr,ylmgq,sfacgq)
! generate the regularised Coulomb Green's function in G+q-space
  call gengclgq(.true.,iq,ngrf,gqc,gclgq)
! use the symmetric form
  gclgq(:)=sqrt(gclgq(:))
! zero the response function (stored in epsi)
  epsi(:,:,:)=0.d0
  do ik=1,nkptnr
! compute v^1/2 chi0 v^1/2
    call genvchi0(.false.,ik,lock,0.d0,vql(:,iq),gclgq,jlgqr,ylmgq,sfacgq, &
     ngrf,epsi)
  end do
! negate and add delta(G,G')
  epsi(:,:,:)=-epsi(:,:,:)
  do ig=1,ngrf
    epsi(ig,ig,:)=epsi(ig,ig,:)+1.d0
  end do
!-------------------------------------!
!     invert epsilon over G-space     !
!-------------------------------------!
  do iw=1,nwrf
    call zminv(ngrf,epsi(:,:,iw))
  end do
! write inverse RPA epsilon to EPSINV.OUT
  call putepsinv(iq,epsi)
! end loop over q-points
end do
deallocate(vgqc,gqc,gclgq,jlgqr)
deallocate(ylmgq,sfacgq,epsi)
return
end subroutine

