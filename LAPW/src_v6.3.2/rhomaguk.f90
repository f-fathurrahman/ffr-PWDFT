
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhomaguk(ik0,lock,evecu)
use modmain
use modulr
use modomp
implicit none
! arguments
integer, intent(in) :: ik0
integer(8), intent(in) :: lock(nqpt)
complex(8), intent(in) :: evecu(nstulr,nstulr)
! local variables
integer ik,ikpa,ist,i,j
integer ngk0,ispn,is,ias
integer npc,ir,nthd
real(8) wm,wi
! automatic arrays
integer idx(nstsv)
! allocatable arrays
integer(8), allocatable :: lockl(:)
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgk(:,:,:)
complex(8), allocatable :: wfrmt(:,:,:,:),wfrgk(:,:,:)
complex(8), allocatable :: wfir(:,:),zfft(:)
! central k-point
ik=(ik0-1)*nkpa+1
! number of G+k-vectors for central k-point
ngk0=ngk(1,ik)
! initialise the local OpenMP locks
allocate(lockl(nqpt))
do ir=1,nqpt
  call omp_init_lock(lockl(ir))
end do
! get the eigenvectors from file
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk0,vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfgk(ngk0,nspinor,nstsv))
call genwfsv(.false.,.true.,nstsv,idx,ngridg,igfft,ngk(1,ik),igkig(:,1,ik), &
 apwalm,evecfv,evecsv,wfmt,ngk0,wfgk)
deallocate(apwalm,evecfv,evecsv)
allocate(wfrmt(npcmtmax,natmtot,nspinor,nqpt),wfrgk(ngk0,nspinor,nqpt))
! loop over ultra long-range states
do j=1,nstulr
  wm=occulr(j,ik0)
  if (abs(wm).lt.epsocc) cycle
  wm=wm*wkpt(ik)
  wi=wm/omega
! zero the ultra long-range wavefunctions
  call holdthd(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP PRIVATE(ir,ispn) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
  do ir=1,nqpt
    do ispn=1,nspinor
      call wfmt0(wfrmt(:,:,ispn,ir))
    end do
  end do
!$OMP SECTION
  wfrgk(:,:,:)=0.d0
!$OMP END PARALLEL SECTIONS
  call freethd(nthd)
! parallel loop over second-variational states
  call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfft,ikpa,i,ir) &
!$OMP NUM_THREADS(nthd)
  allocate(zfft(nqpt))
!$OMP DO
  do ist=1,nstsv
    zfft(:)=0.d0
! loop over kappa-points
    do ikpa=1,nkpa
      i=(ikpa-1)*nstsv+ist
! store the wavefunction in Q-space
      zfft(iqfft(ikpa))=evecu(i,j)
    end do
! Fourier transform to R-space
    call zfftifc(3,ngridq,1,zfft)
! loop over R-points
    do ir=1,nqpt
      call omp_set_lock(lockl(ir))
      call wfadd(zfft(ir),wfmt(:,:,:,ist),wfgk(:,:,ist),wfrmt(:,:,:,ir), &
       wfrgk(:,:,ir))
      call omp_unset_lock(lockl(ir))
    end do
  end do
!$OMP END DO
  deallocate(zfft)
!$OMP END PARALLEL
  call freethd(nthd)
! parallel loop over R-points
  call holdthd(nqpt,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir,ispn,ias,is,npc) &
!$OMP NUM_THREADS(nthd)
  allocate(wfir(ngtot,nspinor))
!$OMP DO
  do ir=1,nqpt
    do ispn=1,nspinor
! Fourier transform the interstitial part to real-space
      wfir(:,ispn)=0.d0
      wfir(igfft(igkig(1:ngk0,1,ik)),ispn)=wfrgk(1:ngk0,ispn,ir)
      call zfftifc(3,ngridg,1,wfir(:,ispn))
    end do
! add to the density and magnetisation
    call omp_set_lock(lock(ir))
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      if (spinpol) then
        if (ncmag) then
          call rmk1(npc,wm,wfrmt(:,ias,1,ir),wfrmt(:,ias,2,ir), &
           rhormt(:,ias,ir),magrmt(:,ias,1,ir),magrmt(:,ias,2,ir), &
           magrmt(:,ias,3,ir))
        else
          call rmk2(npc,wm,wfrmt(:,ias,1,ir),wfrmt(:,ias,2,ir), &
           rhormt(:,ias,ir),magrmt(:,ias,1,ir))
        end if
      else
        call rmk3(npc,wm,wfrmt(:,ias,1,ir),rhormt(:,ias,ir))
      end if
    end do
    if (spinpol) then
      if (ncmag) then
        call rmk1(ngtot,wi,wfir,wfir(:,2),rhorir(:,ir),magrir(:,1,ir), &
         magrir(:,2,ir),magrir(:,3,ir))
      else
        call rmk2(ngtot,wi,wfir,wfir(:,2),rhorir(:,ir),magrir(:,1,ir))
      end if
    else
      call rmk3(ngtot,wi,wfir,rhorir(:,ir))
    end if
    call omp_unset_lock(lock(ir))
! end loop over R-points
  end do
!$OMP END DO
  deallocate(wfir)
!$OMP END PARALLEL
  call freethd(nthd)
! end loop over long-range states
end do
! destroy the local OpenMP locks
do ir=1,nqpt
  call omp_destroy_lock(lockl(ir))
end do
deallocate(lockl)
deallocate(wfmt,wfgk,wfrmt,wfrgk)
return

contains

subroutine wfmt0(wfmt)
implicit none
! arguments
complex(8), intent(out) :: wfmt(npcmtmax,natmtot)
! local variables
integer is,ias
do ias=1,natmtot
  is=idxis(ias)
  wfmt(1:npcmt(is),ias)=0.d0
end do
return
end subroutine

subroutine wfadd(za,wfmt1,wfgk1,wfmt2,wfgk2)
implicit none
! arguments
complex(8), intent(in) :: za
complex(8), intent(in) :: wfmt1(npcmtmax,natmtot,nspinor)
complex(8), intent(in) :: wfgk1(ngk0,nspinor)
complex(8), intent(inout) :: wfmt2(npcmtmax,natmtot,nspinor)
complex(8), intent(inout) :: wfgk2(ngk0,nspinor)
! local variables
integer ispn,is,ias
do ispn=1,nspinor
  do ias=1,natmtot
    is=idxis(ias)
    call zaxpy(npcmt(is),za,wfmt1(:,ias,ispn),1,wfmt2(:,ias,ispn),1)
  end do
end do
do ispn=1,nspinor
  call zaxpy(ngk0,za,wfgk1(:,ispn),1,wfgk2(:,ispn),1)
end do
return
end subroutine

pure subroutine rmk1(n,wo,wf1,wf2,rho,mag1,mag2,mag3)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag1(n),mag2(n),mag3(n)
! local variables
integer i
real(8) wo2,t1,t2
complex(8) z1,z2
wo2=2.d0*wo
do i=1,n
  z1=wf1(i)
  z2=wf2(i)
  t1=dble(z1)**2+aimag(z1)**2
  t2=dble(z2)**2+aimag(z2)**2
  z1=conjg(z1)*z2
  rho(i)=rho(i)+wo*(t1+t2)
  mag1(i)=mag1(i)+wo2*dble(z1)
  mag2(i)=mag2(i)+wo2*aimag(z1)
  mag3(i)=mag3(i)+wo*(t1-t2)
end do
return
end subroutine

pure subroutine rmk2(n,wo,wf1,wf2,rho,mag)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag(n)
! local variables
integer i
real(8) t1,t2
do i=1,n
  t1=dble(wf1(i))**2+aimag(wf1(i))**2
  t2=dble(wf2(i))**2+aimag(wf2(i))**2
  rho(i)=rho(i)+wo*(t1+t2)
  mag(i)=mag(i)+wo*(t1-t2)
end do
return
end subroutine

pure subroutine rmk3(n,wo,wf,rho)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf(n)
real(8), intent(inout) :: rho(n)
rho(:)=rho(:)+wo*(dble(wf(:))**2+aimag(wf(:))**2)
return
end subroutine

end subroutine

