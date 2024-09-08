
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potksu
use modmain
use modulr
use modomp
implicit none
! local variables
integer ifq,idm,is,ias
integer npc,nrc,nrci,nthd
! allocatable arrays
real(8), allocatable :: rfmt1(:),rfmt2(:)
! compute the ultra long-range Coulomb potential
call potcoulu
! compute the ultra long-range exchange-correlation potential and fields
call potxcu
! subtract the normal Kohn-Sham potential for Q=0
call holdthd(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt1,rfmt2,is) &
!$OMP PRIVATE(nrc,nrci,npc) &
!$OMP NUM_THREADS(nthd)
allocate(rfmt1(npcmtmax),rfmt2(npcmtmax))
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  call rfmtftoc(nrc,nrci,vsmt(:,ias),rfmt1)
  call rbsht(nrc,nrci,rfmt1,rfmt2)
  vsqmt(1:npc,ias,1)=vsqmt(1:npc,ias,1)-rfmt2(1:npc)
end do
!$OMP END DO
deallocate(rfmt1,rfmt2)
!$OMP END PARALLEL
call freethd(nthd)
vsqir(:,1)=vsqir(:,1)-vsir(:)
! multiply vsqir by the characteristic function
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ifq=1,nfqrz
  vsqir(:,ifq)=vsqir(:,ifq)*cfunir(:)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
if (.not.spinpol) return
! subtract the normal Kohn-Sham magnetic field for Q=0 in the muffin-tins
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    npc=npcmt(is)
    bsqmt(1:npc,ias,idm,1)=bsqmt(1:npc,ias,idm,1)-bsmt(1:npc,ias,idm)
  end do
end do
! multiply bsqir by the characteristic function
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(idm) &
!$OMP NUM_THREADS(nthd)
do ifq=1,nfqrz
  do idm=1,ndmag
    bsqir(:,idm,ifq)=bsqir(:,idm,ifq)*cfunir(:)
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! subtract the normal Kohn-Sham magnetic field for Q=0 in the interstitial
! (this is already multiplied by the characteristic function)
do idm=1,ndmag
  bsqir(:,idm,1)=bsqir(:,idm,1)-bsir(:,idm)
end do
! reduce the external magnetic field if required
if (reducebf.lt.1.d0) then
  bfcq(:,:)=bfcq(:,:)*reducebf
  bfcmtq(:,:,:)=bfcmtq(:,:,:)*reducebf
end if
return
end subroutine

