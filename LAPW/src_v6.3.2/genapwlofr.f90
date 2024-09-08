
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genapwlofr
use modomp
implicit none
! local variables
integer nthd
call holdthd(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
! generate the APW radial functions
call genapwfr
!$OMP SECTION
! generate the local-orbital radial functions
call genlofr
!$OMP END PARALLEL SECTIONS
call freethd(nthd)
call holdthd(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
! compute the overlap radial integrals
call olprad
!$OMP SECTION
! compute the Hamiltonian radial integrals
call hmlrad
!$OMP END PARALLEL SECTIONS
call freethd(nthd)
return
end subroutine

