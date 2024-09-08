
! Copyright (C) 2007-2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genkmat
! !INTERFACE:
subroutine genkmat(tfv,tvclcr)
! !USES:
use modmain
use modmpi
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   tfv    : .true. if the matrix elements are to be expressed in the
!            first-variational basis; second-variational otherwise (in,logical)
!   tvclvr : .true. if the non-local Coulomb potential from the core states is
!            to be included in the kinetic matrix elements (in,logical)
! !DESCRIPTION:
!   Computes the kinetic matrix elements in the first- or second-variational
!   basis and stores them in the file {\tt KMAT.OUT}. See routine {\tt putkmat}.
!
! !REVISION HISTORY:
!   Created January 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tfv,tvclcr
! local variables
integer ik,is,ias
integer nrc,nrci,nthd
! allocatable arrays
real(8), allocatable :: vmt(:,:),vir(:),rfmt(:)
allocate(vmt(npcmtmax,natmtot),vir(ngtot))
! convert muffin-tin Kohn-Sham potential to spherical coordinates on a coarse
! radial mesh
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
! multiply Kohn-Sham interstitial potential by characteristic function
vir(:)=vsir(:)*cfunir(:)
if (mp_mpi) write(*,*)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! loop over k-points
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL(genkmat_)
  write(*,'("Info(genkmat): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL(genkmat_)
  call putkmat(tfv,tvclcr,ik,vmt,vir)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
deallocate(vmt,vir)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
return
end subroutine
!EOC

