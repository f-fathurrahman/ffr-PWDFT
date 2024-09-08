
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genws
use modmain
use modomp
implicit none
! local variables
integer is,ias,nthd
integer nrc,nrci
! allocatable arrays
real(8), allocatable :: rfmt(:)
if (xcgrad.ne.4) return
! muffin-tin effective tau-DFT potential
call holdthd(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nrc,nrci) &
!$OMP NUM_THREADS(nthd)
allocate(rfmt(npcmtmax))
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
! convert to coarse radial mesh and spherical coordinates
  call rfmtftoc(nrc,nrci,wxcmt(:,ias),rfmt)
  call rbsht(nrc,nrci,rfmt,wsmt(:,ias))
end do
!$OMP END DO
deallocate(rfmt)
!$OMP END PARALLEL
call freethd(nthd)
! interstitial tau-DFT potential
wsir(:)=wxcir(:)*cfunir(:)
return
end subroutine

