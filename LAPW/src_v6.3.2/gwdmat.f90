
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gwdmat
use modmain
use modgw
use modmpi
use modomp
implicit none
! local variables
integer ik,nthd
! initialise universal variables
call init0
call init1
call init3
! read Fermi energy from file
call readfermi
! get the eigenvalues from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
end do
! determine the GW Fermi energy
call gwefermi
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
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(gwdmat):")')
  write(*,'(" GW density matrices determined for each k-point")')
  write(*,*)
  write(*,'(" Natural orbitals and occupation numbers written to")')
  write(*,'(" EVECSV.OUT and OCCSV.OUT, respectively")')
end if
return
end subroutine

