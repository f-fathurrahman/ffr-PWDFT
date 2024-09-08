
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine timestep
use modmain
use modtddft
use modmpi
use modomp
implicit none
! local variables
integer ik,is,ias,i,j
integer nrc,nrci,nthd
real(8) dt,t1
complex(8) z1
! allocatable arrays
real(8), allocatable :: vmt(:,:),vir(:),rfmt(:),w(:)
complex(8), allocatable :: evecsv(:,:),evectv(:,:),evecsvt(:,:)
complex(8), allocatable :: a(:,:),b(:,:),c(:,:)
if (itimes.ge.ntimes) then
  write(*,*)
  write(*,'("Error(timestep): itimes >= ntimes : ",2I8)') itimes,ntimes
  write(*,*)
  stop
end if
! convert muffin-tin Kohn-Sham potential to spherical coordinates
allocate(vmt(npcmtmax,natmtot))
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
  call rfmtftoc(nrc,nrci,vsmt(:,ias),rfmt)
  call rbsht(nrc,nrci,rfmt,vmt(:,ias))
end do
!$OMP END DO
deallocate(rfmt)
!$OMP END PARALLEL
call freethd(nthd)
! multiply interstitial potential by characteristic function
allocate(vir(ngtot))
vir(:)=vsir(:)*cfunir(:)
! backup existing time-dependent Kohn-Sham eigenvectors if required
call tdbackup
! loop over k-points
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecsv,evectv,evecsvt) &
!$OMP PRIVATE(w,a,b,c,i,j,dt,t1,z1) &
!$OMP NUM_THREADS(nthd)
allocate(evecsv(nstsv,nstsv),evectv(nstsv,nstsv),evecsvt(nstsv,nstsv))
allocate(w(nstsv),a(nstsv,nstsv),b(nstsv,nstsv),c(nstsv,nstsv))
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! generate the Hamiltonian matrix in the ground-state second-variational basis
  call genhmlt(ik,vmt,vir,evectv)
! diagonalise the Hamiltonian to get third-variational eigenvectors
  if (spinpol.and.(.not.ncmag)) then
! collinear case requires block diagonalisation
    call eveqnz(nstfv,nstsv,evectv,w)
    i=nstfv+1
    call eveqnz(nstfv,nstsv,evectv(i,i),w(i))
    do i=1,nstfv
      do j=1,nstfv
        evectv(i,j+nstfv)=0.d0
        evectv(i+nstfv,j)=0.d0
      end do
    end do
  else
! non-collinear or spin-unpolarised: full diagonalisation
    call eveqnz(nstsv,nstsv,evectv,w)
  end if
! read in ground-state eigenvectors
  call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! convert third-variational eigenvectors to first-variational basis
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,evectv,nstsv,zzero,a, &
   nstsv)
! time propagate instantaneous eigenvectors across one time step
  dt=times(itimes+1)-times(itimes)
  if (tdphi.eq.0.d0) then
! real time evolution
    do i=1,nstsv
      t1=-w(i)*dt
      z1=cmplx(cos(t1),sin(t1),8)
      b(:,i)=z1*a(:,i)
    end do
  else
! complex time evolution
    do i=1,nstsv
      t1=-w(i)*dt
      z1=t1*cmplx(sin(tdphi),cos(tdphi),8)
      z1=exp(z1)
      b(:,i)=z1*a(:,i)
    end do
  end if
! read in time-dependent Kohn-Sham eigenvectors (first-variational basis)
  call getevecsv(filext,ik,vkl(:,ik),evecsvt)
! apply time-evolution operator
  call zgemm('C','N',nstsv,nstsv,nstsv,zone,a,nstsv,evecsvt,nstsv,zzero,c,nstsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,b,nstsv,c,nstsv,zzero,evecsvt,nstsv)
! orthonormalise the eigenvectors if required
  if (tdphi.ne.0.d0) call orthevsv(evecsvt)
! write the new eigenvectors to file
  call putevecsv(filext,ik,evecsvt)
end do
!$OMP END DO
deallocate(evecsv,evectv,evecsvt)
deallocate(w,a,b,c)
!$OMP END PARALLEL
call freethd(nthd)
deallocate(vmt,vir)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
return
end subroutine

