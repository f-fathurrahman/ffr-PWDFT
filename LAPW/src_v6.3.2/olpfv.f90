
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpfv(nmatp,ngp,igpig,apwalm,o)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: nmatp,ngp,igpig(ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: o(nmatp,nmatp)
! local variables
integer ias,i
integer nthd1,nthd2
! zero the upper triangular part of the matrix
do i=1,nmatp
  o(1:i,i)=0.d0
end do
call holdthd(2,nthd1)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP PRIVATE(ias) &
!$OMP NUM_THREADS(nthd1)
!$OMP SECTION
do ias=1,natmtot
  call olpaa(tefvr,ias,ngp,apwalm(:,:,:,ias),nmatp,o)
end do
call olpistl(ngp,igpig,nmatp,o)
!$OMP SECTION
call holdthd(natmtot,nthd2)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd2)
do ias=1,natmtot
  call olpalo(ias,ngp,apwalm(:,:,:,ias),nmatp,o)
  call olplolo(ias,ngp,nmatp,o)
end do
!$OMP END PARALLEL DO
call freethd(nthd2)
!$OMP END PARALLEL SECTIONS
call freethd(nthd1)
return
end subroutine

