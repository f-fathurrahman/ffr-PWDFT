
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dpotks
use modmain
use modphonon
use modomp
implicit none
! local variables
integer is,ias,nthd
integer nr,nri,np
! allocatable arrays
complex(8), allocatable :: zfmt(:)
! convert density derivative to spherical coordinates
call holdthd(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt,is,nr,nri,np) &
!$OMP NUM_THREADS(nthd)
allocate(zfmt(npmtmax))
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  zfmt(1:np)=drhomt(1:np,ias)
  call zbsht(nr,nri,zfmt,drhomt(:,ias))
end do
!$OMP END DO
deallocate(zfmt)
!$OMP END PARALLEL
call freethd(nthd)
! compute the exchange-correlation potential derivative
call dpotxc
! convert density derivative to spherical harmonics
call holdthd(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt,is,nr,nri,np) &
!$OMP NUM_THREADS(nthd)
allocate(zfmt(npmtmax))
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  zfmt(1:np)=drhomt(1:np,ias)
  call zfsht(nr,nri,zfmt,drhomt(:,ias))
end do
!$OMP END DO
deallocate(zfmt)
!$OMP END PARALLEL
call freethd(nthd)
! generate the Coulomb potential derivative
call dpotcoul
! add to the Kohn-Sham potential derivative
do ias=1,natmtot
  is=idxis(ias)
  np=npmt(is)
  dvsmt(1:np,ias)=dvsmt(1:np,ias)+dvclmt(1:np,ias)
end do
call zaxpy(ngtot,zone,dvclir,1,dvsir,1)
! remove the gradient part of the potential derivative for displaced muffin-tin
np=npmt(isph)
dvsmt(1:np,iasph)=dvsmt(1:np,iasph)+gvsmt(1:np)
return
end subroutine

