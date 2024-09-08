
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potcoul
! !INTERFACE:
subroutine potcoul
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Calculates the Coulomb potential of the real charge density stored in the
!   global variables {\tt rhomt} and {\tt rhoir} by solving Poisson's equation.
!   These variables are coverted to complex representations and passed to the
!   routine {\tt zpotcoul}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ias,nthd
integer nr,nri,ir,i
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:),zvclir(:)
allocate(zrhomt(npmtmax,natmtot))
! convert real muffin-tin charge density to complex spherical harmonic expansion
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call rtozfmt(nrmt(is),nrmti(is),rhomt(:,ias),zrhomt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! solve the complex Poisson's equation in the muffin-tins
allocate(zvclmt(npmtmax,natmtot))
call genzvclmt(nrmt,nrmti,nrmtmax,rlmt,wprmt,npmtmax,zrhomt,zvclmt)
deallocate(zrhomt)
! add the nuclear monopole potentials
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  i=1
  do ir=1,nri
    zvclmt(i,ias)=zvclmt(i,ias)+vcln(ir,is)
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    zvclmt(i,ias)=zvclmt(i,ias)+vcln(ir,is)
    i=i+lmmaxo
  end do
end do
! store real interstitial charge density in complex array
allocate(zrhoir(ngtot))
zrhoir(:)=rhoir(:)
! solve Poisson's equation in the entire unit cell
allocate(zvclir(ngtot))
call zpotcoul(nrmt,nrmti,npmt,npmti,nrmtmax,rlmt,ngridg,igfft,ngvec,gc,gclg, &
 ngvec,jlgrmt,ylmg,sfacg,zrhoir,npmtmax,zvclmt,zvclir)
! convert complex muffin-tin potential to real spherical harmonic expansion
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  call ztorfmt(nrmt(is),nrmti(is),zvclmt(:,ias),vclmt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! store complex interstitial potential in real array
vclir(:)=dble(zvclir(:))
deallocate(zrhoir,zvclmt,zvclir)
! apply constant electric field if required
if (tefield) call potefield
return
end subroutine
!EOC

