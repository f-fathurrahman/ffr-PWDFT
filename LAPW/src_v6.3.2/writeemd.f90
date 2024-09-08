
! Copyright (C) 2012 S. Dugdale, D. Ernsting and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeemd
use modmain
use modpw
use modmpi
use modomp
implicit none
! local variables
integer ik,ihk,recl
integer ist,ispn,nthd
real(8) sum,t1
complex(8) z1
! allocatable arrays
real(8), allocatable :: emd(:)
complex(8), allocatable :: wfpw(:,:,:)
if (spinsprl) then
  write(*,*)
  write(*,'("Error(writeemd): electron momentum density not available for &
   &spin-spirals")')
  write(*,*)
  stop
end if
! initialise universal variables
call init0
call init1
call init4
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the occupancies from file
do ik=1,nkpt
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
! delete existing EMD.OUT
if (mp_mpi) then
  open(160,file='EMD.OUT')
  close(160,status='DELETE')
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
allocate(emd(nhkmax))
inquire(iolength=recl) vkl(:,1),nhkmax,emd
deallocate(emd)
open(160,file='EMD.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
! loop over k-points
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(emd,wfpw,ihk,sum) &
!$OMP PRIVATE(ist,ispn,z1,t1) &
!$OMP NUM_THREADS(nthd)
allocate(emd(nhkmax),wfpw(nhkmax,nspinor,nstsv))
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(writeemd_)
  write(*,'("Info(writeemd): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(writeemd_)
! Fourier transform the wavefunctions
  call genwfpw(vkl(:,ik),ngk(1,ik),igkig(:,1,ik),vgkl(:,:,1,ik), &
   vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),nhk(1,ik),vhkc(:,:,1,ik), &
   hkc(:,1,ik),sfachk(:,:,1,ik),wfpw)
! loop over all H+k-vectors
  do ihk=1,nhk(1,ik)
! sum over occupied states and spins
    sum=0.d0
    do ist=1,nstsv
      do ispn=1,nspinor
        z1=wfpw(ihk,ispn,ist)
        t1=dble(z1)**2+aimag(z1)**2
        sum=sum+occsv(ist,ik)*t1
      end do
    end do
    emd(ihk)=sum
  end do
!$OMP CRITICAL(u160)
  write(160,rec=ik) vkl(:,ik),nhk(1,ik),emd
!$OMP END CRITICAL(u160)
end do
!$OMP END DO
deallocate(emd,wfpw)
!$OMP END PARALLEL
call freethd(nthd)
close(160)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(writeemd): electron momentum density written to EMD.OUT")')
  write(*,'(" for all H+k-vectors up to |H+k| < hkmax")')
end if
return
end subroutine

