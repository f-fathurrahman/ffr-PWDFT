
! Copyright (C) 2018 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gwrhomag
use modmain
use modgw
use modmpi
use modomp
implicit none
! local variables
integer ik,nthd
! compute the GW density matrices and write the natural orbitals and occupation
! numbers to EVECSV.OUT and OCCSV.OUT, respectively
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  call gwdmatk(ik)
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! determine the density and magnetisation in the usual way
call rhomag
return
end subroutine

