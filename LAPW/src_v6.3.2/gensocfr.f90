
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine gensocfr
use modmain
use modomp
implicit none
! local variables
integer is,ias,nthd
integer nr,nri,ir,irc
real(8) cso,rm
! allocatable arrays
real(8), allocatable :: vr(:),dvr(:)
if (.not.spinorb) return
! coefficient of spin-orbit coupling
cso=socscf/(4.d0*solsc**2)
call holdthd(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vr,dvr,is,nr,nri) &
!$OMP PRIVATE(irc,ir,rm) &
!$OMP NUM_THREADS(nthd)
allocate(vr(nrmtmax),dvr(nrmtmax))
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
! radial derivative of the spherical part of the Kohn-Sham potential
  call rfmtlm(1,nr,nri,vsmt(:,ias),vr)
  vr(1:nr)=vr(1:nr)*y00
  call fderiv(1,nr,rlmt(:,1,is),vr,dvr)
  irc=0
  do ir=1,nr,lradstp
    irc=irc+1
    rm=1.d0-2.d0*cso*vr(ir)
    socfr(irc,ias)=cso*dvr(ir)/(rsp(ir,is)*rm**2)
  end do
end do
!$OMP END DO
deallocate(vr,dvr)
!$OMP END PARALLEL
call freethd(nthd)
return
end subroutine

