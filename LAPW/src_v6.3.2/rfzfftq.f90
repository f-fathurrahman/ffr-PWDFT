
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfzfftq(sgn,nf,rfmt,rfir,zfmt,zfir)
use modmain
use modmpi
use modomp
implicit none
! arguments
integer, intent(in) :: sgn,nf
real(8), intent(inout) :: rfmt(npcmtmax,natmtot,nf,nqpt)
real(8), intent(inout) :: rfir(ngtot,nf,nqpt)
complex(8), intent(inout) :: zfmt(npcmtmax,natmtot,nf,nfqrz)
complex(8), intent(inout) :: zfir(ngtot,nf,nfqrz)
! local variables
integer jf,is,ias,ir
integer npc,i,n,nthd
if (sgn.eq.-1) then
! loop over the number of functions
  do jf=1,nf
! Fourier transform the muffin-tin function
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      if (np_mpi.gt.1) zfmt(1:npc,ias,jf,:)=0.d0
      call holdthd(npc/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
      do i=1,npc
! distribute among MPI processes
        if (mod(i-1,np_mpi).ne.lp_mpi) cycle
        call rzfftifc(3,ngridq,-1,rfmt(i,ias,jf,:),zfmt(i,ias,jf,:))
      end do
!$OMP END PARALLEL DO
      call freethd(nthd)
    end do
! Fourier transform the interstitial function
    if (np_mpi.gt.1) zfir(:,jf,:)=0.d0
    call holdthd(ngtot/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
    do ir=1,ngtot
! distribute among MPI processes
      if (mod(ir-1,np_mpi).ne.lp_mpi) cycle
      call rzfftifc(3,ngridq,-1,rfir(ir,jf,:),zfir(ir,jf,:))
    end do
!$OMP END PARALLEL DO
    call freethd(nthd)
! end loop over number of functions
  end do
  if (np_mpi.gt.1) then
    n=npcmtmax*natmtot*nf*nfqrz
    call mpi_allreduce(mpi_in_place,zfmt,n,mpi_double_complex,mpi_sum,mpicom, &
     ierror)
    n=ngtot*nf*nfqrz
    call mpi_allreduce(mpi_in_place,zfir,n,mpi_double_complex,mpi_sum,mpicom, &
     ierror)
  end if
else
! loop over the number of functions
  do jf=1,nf
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      if (np_mpi.gt.1) rfmt(1:npc,ias,jf,:)=0.d0
      call holdthd(npc/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
      do i=1,npc
! distribute among MPI processes
        if (mod(i-1,np_mpi).ne.lp_mpi) cycle
        call rzfftifc(3,ngridq,1,rfmt(i,ias,jf,:),zfmt(i,ias,jf,:))
      end do
!$OMP END PARALLEL DO
      call freethd(nthd)
    end do
    if (np_mpi.gt.1) rfir(:,jf,:)=0.d0
    call holdthd(ngtot/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
    do ir=1,ngtot
! distribute among MPI processes
      if (mod(ir-1,np_mpi).ne.lp_mpi) cycle
      call rzfftifc(3,ngridq,1,rfir(ir,jf,:),zfir(ir,jf,:))
    end do
!$OMP END PARALLEL DO
    call freethd(nthd)
! end loop over number of functions
  end do
  if (np_mpi.gt.1) then
    n=npcmtmax*natmtot*nf*nqpt
    call mpi_allreduce(mpi_in_place,rfmt,n,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
    n=ngtot*nf*nqpt
    call mpi_allreduce(mpi_in_place,rfir,n,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
  end if
end if
return
end subroutine

