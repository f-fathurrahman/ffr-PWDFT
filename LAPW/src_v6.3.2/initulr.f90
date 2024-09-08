
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initulr
use modmain
use modulr
use modomp
implicit none
! local variables
integer ik0,ik,ist,jst
integer iq,ifq,ig,i,nthd
real(8) t1
! allocatable arrays
integer, allocatable :: idx(:)
! allocate long-range density and magnetisation arrays
if (allocated(rhormt)) deallocate(rhormt)
allocate(rhormt(npcmtmax,natmtot,nqpt))
if (allocated(rhorir)) deallocate(rhorir)
allocate(rhorir(ngtot,nqpt))
if (allocated(magrmt)) deallocate(magrmt)
if (allocated(magrir)) deallocate(magrir)
if (spinpol) then
  allocate(magrmt(npcmtmax,natmtot,ndmag,nqpt))
  allocate(magrir(ngtot,ndmag,nqpt))
end if
if (allocated(rhoqmt)) deallocate(rhoqmt)
allocate(rhoqmt(npcmtmax,natmtot,nfqrz))
if (allocated(rhoqir)) deallocate(rhoqir)
allocate(rhoqir(ngtot,nfqrz))
if (allocated(chgmtru)) deallocate(chgmtru)
allocate(chgmtru(natmtot,nqpt))
if (allocated(magqmt)) deallocate(magqmt)
if (allocated(magqir)) deallocate(magqir)
if (allocated(mommtru)) deallocate(mommtru)
if (spinpol) then
  allocate(magqmt(npcmtmax,natmtot,ndmag,nfqrz))
  allocate(magqir(ngtot,ndmag,nfqrz))
  allocate(mommtru(ndmag,natmtot,nqpt))
end if
! allocate potential and magnetic field arrays
if (allocated(vclq)) deallocate(vclq)
allocate(vclq(nfqrz))
if (allocated(vsqmt)) deallocate(vsqmt)
allocate(vsqmt(npcmtmax,natmtot,nfqrz))
if (allocated(vsqir)) deallocate(vsqir)
allocate(vsqir(ngtot,nfqrz))
if (allocated(bfcq)) deallocate(bfcq)
if (allocated(bfcmtq)) deallocate(bfcmtq)
if (allocated(bsqmt)) deallocate(bsqmt)
if (allocated(bsqir)) deallocate(bsqir)
if (spinpol) then
  allocate(bfcq(ndmag,nfqrz))
  allocate(bfcmtq(natmtot,ndmag,nfqrz))
  allocate(bsqmt(npcmtmax,natmtot,ndmag,nfqrz))
  allocate(bsqir(ngtot,ndmag,nfqrz))
end if
! G+Q-vector arrays
if (allocated(vgqc)) deallocate(vgqc)
allocate(vgqc(3,ngvec,nfqrz))
if (allocated(gqc)) deallocate(gqc)
allocate(gqc(ngvec,nfqrz))
if (allocated(ylmgq)) deallocate(ylmgq)
allocate(ylmgq(lmmaxo,ngvec,nfqrz))
if (allocated(sfacgq)) deallocate(sfacgq)
allocate(sfacgq(ngvec,natmtot,nfqrz))
if (allocated(gclq)) deallocate(gclq)
allocate(gclq(nqpt))
if (allocated(gclgq)) deallocate(gclgq)
allocate(gclgq(ngvec,nfqrz))
if (allocated(jlgqrmt)) deallocate(jlgqrmt)
allocate(jlgqrmt(0:lnpsd,ngvec,nspecies,nfqrz))
! find the maximum size of the spherical Bessel function array over all species
call findnjcmax
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(iq,ig,t1) &
!$OMP NUM_THREADS(nthd)
do ifq=1,nfqrz
  iq=iqrzf(ifq)
  do ig=1,ngvec
! determine the G+Q-vectors
    vgqc(:,ig,ifq)=vgc(:,ig)+vqc(:,iq)
! G+Q-vector length
    gqc(ig,ifq)=sqrt(vgqc(1,ig,ifq)**2+vgqc(2,ig,ifq)**2+vgqc(3,ig,ifq)**2)
! spherical harmonics for G+Q-vectors
    call genylmv(lmaxo,vgqc(:,ig,ifq),ylmgq(:,ig,ifq))
  end do
! structure factors for G+Q-vectors
  call gensfacgp(ngvec,vgqc(:,:,ifq),ngvec,sfacgq(:,:,ifq))
! generate the Coulomb Green's function in Q-space with small Q cut-off
  t1=sqrt(vqc(1,iq)**2+vqc(2,iq)**2+vqc(3,iq)**2)
  if (t1.gt.q0cut+epslat) then
    gclq(iq)=fourpi/t1**2
  else
    gclq(iq)=0.d0
  end if
! generate the Coulomb Green's function in G+Q-space
  call gengclgq(.true.,iq,ngvec,gqc(:,ifq),gclgq(:,ifq))
! compute the spherical Bessel functions j_l(|G+Q|R_mt)
  call genjlgprmt(lnpsd,ngvec,gqc(:,ifq),ngvec,jlgqrmt(:,:,:,ifq))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! number of long-range states
nstulr=nstsv*nkpa
! allocate eigenvalue array
if (allocated(evalu)) deallocate(evalu)
allocate(evalu(nstulr,nkpt0))
! allocate the occupation number array
if (allocated(occulr)) deallocate(occulr)
allocate(occulr(nstulr,nkpt0))
! initialise the occupation numbers
allocate(idx(nstulr))
do ik0=1,nkpt0
  ik=(ik0-1)*nkpa+1
  call sortidx(nstulr,occsv(1,ik),idx)
  do ist=1,nstulr
    i=idx(nstulr-ist+1)-1
    ik=(ik0-1)*nkpa+i/nstsv+1
    jst=mod(i,nstsv)+1
    occulr(ist,ik0)=occsv(jst,ik)
  end do
end do
deallocate(idx)
return
end subroutine

