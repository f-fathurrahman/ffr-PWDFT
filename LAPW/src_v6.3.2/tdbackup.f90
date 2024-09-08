
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tdbackup
use modmain
use modtddft
use modmpi
implicit none
! local variables
integer ik
! allocatable arrays
complex(8), allocatable :: evecsvt(:,:)
if (ntsbackup.le.0) return
if (mod(itimes-1,ntsbackup).ne.0) return
if (mp_mpi) then
  allocate(evecsvt(nstsv,nstsv))
  do ik=1,nkpt
! read in time-dependent Kohn-Sham eigenvectors
    call getevecsv(filext,ik,vkl(:,ik),evecsvt)
! write eigenvectors to backup file
    call putevecsv('_TD_BACKUP.OUT',ik,evecsvt)
  end do
  deallocate(evecsvt)
! write the time-step backup file
  open(50,file='TIMESTEP_BACKUP.OUT',form='FORMATTED')
  write(50,'(I8,G18.10)') itimes,times(itimes)
  close(50)
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
return
end subroutine

