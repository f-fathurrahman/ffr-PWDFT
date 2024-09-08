
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genvclijji
! !INTERFACE:
subroutine genvclijji(ikp,vclijji)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced set (in,integer)
!   vclijji : Coulomb matrix elements (out,real(nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Calculates the Coulomb matrix elements of the type $(i-jj-i)$.
!
! !REVISION HISTORY:
!   Created June 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(out) :: vclijji(nstsv,nstsv,nkpt)
! local variables
integer ik,ist1,ist2
integer iv(3),iq,ig,i
real(8) vc(3)
complex(8) z1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqrmt(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(8), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
! external functions
complex(8) zfinp
external zfinp
! allocate local arrays
allocate(vgqc(3,ngvc),gqc(ngvc),gclgq(ngvc))
allocate(jlgqrmt(0:lnpsd,ngvc,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxo,ngvc),sfacgq(ngvc,natmtot))
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv),wfir1(ngtc,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv),wfir2(ngtc,nspinor,nstsv))
allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtc))
allocate(zvclmt(npcmtmax,natmtot),zvclir(ngtc))
! get the eigenvectors from file for non-reduced k-point ikp
call getevecfv(filext,0,vkl(:,ikp),vgkl(:,:,1,ikp),evecfv)
call getevecsv(filext,0,vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),vgkc(:,:,1,ikp),gkc(:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! index to all states
do ist1=1,nstsv
  idx(ist1)=ist1
end do
! calculate the wavefunctions for all states of passed non-reduced k-point ikp
call genwfsv(.false.,.false.,nstsv,idx,ngdc,igfc,ngk(1,ikp),igkig(:,1,ikp), &
 apwalm,evecfv,evecsv,wfmt2,ngtc,wfir2)
! start loop over reduced k-point set
do ik=1,nkpt
! determine the q-vector
  iv(:)=ivk(:,ik)-ivk(:,ikp)
  iv(:)=modulo(iv(:),ngridk(:))
! check if the q-point is in user-defined set
  iv(:)=iv(:)*ngridq(:)
  do i=1,3
    if (modulo(iv(i),ngridk(i)).ne.0) goto 10
  end do
  iv(:)=iv(:)/ngridk(:)
  iq=iqmap(iv(1),iv(2),iv(3))
  vc(:)=vkc(:,ik)-vkc(:,ikp)
  do ig=1,ngvc
! determine the G+q-vectors
    vgqc(:,ig)=vgc(:,ig)+vc(:)
! G+q-vector length
    gqc(ig)=sqrt(vgqc(1,ig)**2+vgqc(2,ig)**2+vgqc(3,ig)**2)
! spherical harmonics for G+q-vectors
    call genylmv(lmaxo,vgqc(:,ig),ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvc,vgqc,ngvc,sfacgq)
! generate the regularised Coulomb Green's function in G+q-space
  call gengclgq(.true.,iq,ngvc,gqc,gclgq)
! compute the required spherical Bessel functions
  call genjlgprmt(lnpsd,ngvc,gqc,ngvc,jlgqrmt)
! find the matching coefficients
  call match(ngk(1,ik),vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
! get the eigenvectors from file
  call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
! calculate the wavefunctions for all states of the reduced k-point
  call genwfsv(.false.,.false.,nstsv,idx,ngdc,igfc,ngk(1,ik),igkig(:,1,ik), &
   apwalm,evecfv,evecsv,wfmt1,ngtc,wfir1)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
  do ist1=1,nstsv
    do ist2=1,nstsv
! calculate the complex overlap density
      call genzrho(.true.,.true.,ngtc,wfmt2(:,:,:,ist2),wfir2(:,:,ist2), &
       wfmt1(:,:,:,ist1),wfir1(:,:,ist1),zrhomt,zrhoir)
! compute the potential and G=0 coefficient of the density
      call genzvclmt(nrcmt,nrcmti,nrcmtmax,rlcmt,wprcmt,npcmtmax,zrhomt,zvclmt)
      call zpotcoul(nrcmt,nrcmti,npcmt,npcmti,nrcmtmax,rlcmt,ngdc,igfc,ngvc, &
       gqc,gclgq,ngvc,jlgqrmt,ylmgq,sfacgq,zrhoir,npcmtmax,zvclmt,zvclir)
      z1=zfinp(zrhomt,zrhoir,zvclmt,zvclir)
      vclijji(ist1,ist2,ik)=wqptnr*dble(z1)
! end loop over ist2
    end do
! end loop over ist1
  end do
10 continue
! end loop over reduced k-point set
end do
deallocate(vgqc,gqc,gclgq,jlgqrmt)
deallocate(apwalm,evecfv,evecsv,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zrhomt,zrhoir,zvclmt,zvclir)
return
end subroutine
!EOC

