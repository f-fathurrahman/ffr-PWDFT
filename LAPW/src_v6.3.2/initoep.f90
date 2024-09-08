
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine initoep
use modmain
implicit none
! local variables
integer ist,ic,m
integer is,ias,nrc,nrci
integer xctype_(3)
! allocatable arrays
real(8), allocatable :: rfmt(:)
! external functions
real(8) rfint
external rfint
! find maximum core states over all species
ncrmax=0
do is=1,nspecies
  ic=0
  do ist=1,nstsp(is)
    if (spcore(ist,is)) then
      do m=-ksp(ist,is),ksp(ist,is)-1
        ic=ic+1
      end do
    end if
  end do
  ncrmax=max(ncrmax,ic)
end do
! allocate the exchange potential and magnetic field
if (allocated(vxmt)) deallocate(vxmt)
allocate(vxmt(npcmtmax,natmtot))
if (allocated(vxir)) deallocate(vxir)
allocate(vxir(ngtot))
if (spinpol) then
  if (allocated(bxmt)) deallocate(bxmt)
  allocate(bxmt(npcmtmax,natmtot,ndmag))
  if (allocated(bxir)) deallocate(bxir)
  allocate(bxir(ngtot,ndmag))
end if
! initialise the exchange potential to LDA
xctype_(1)=3
xctype_(2:)=0
call potxc(.true.,xctype_,rhomt,rhoir,magmt,magir,taumt,tauir,exmt,exir,ecmt, &
 ecir,vxcmt,vxcir,bxcmt,bxcir,wxcmt,wxcir)
allocate(rfmt(npcmtmax))
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  call rfmtftoc(nrc,nrci,vxcmt(:,ias),rfmt)
  call rbsht(nrc,nrci,rfmt,vxmt(:,ias))
end do
deallocate(rfmt)
vxir(:)=vxcir(:)
! determine the constant part of V_xc
vxc0=rfint(vxcmt,vxcir)/omega
return
end subroutine

