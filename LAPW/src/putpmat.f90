
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putpmat(ik)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
! local variables
integer ist,ispn,recl,i
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgk(:,:,:),pmat(:,:,:)
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
! get the eigenvectors from file
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
do ispn=1,nspnfv
  call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
! calculate the wavefunctions for all states
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfgk(ngkmax,nspinor,nstsv))
call genwfsv(.true.,.true.,nstsv,idx,ngridg,igfft,ngk(:,ik),igkig(:,:,ik), &
 apwalm,evecfv,evecsv,wfmt,ngkmax,wfgk)
deallocate(evecfv,evecsv,apwalm)
! calculate the momentum matrix elements
allocate(pmat(nstsv,nstsv,3))
call genpmatk(ngk(:,ik),igkig(:,:,ik),vgkc(:,:,:,ik),wfmt,wfgk,pmat)
deallocate(wfmt,wfgk)
! determine the record length
inquire(iolength=recl) vkl(:,1),nstsv,pmat
! write the matrix elements in the second-variational basis
do i=1,2
  open(150,file='PMAT.OUT',form='UNFORMATTED',access='DIRECT',recl=recl,err=10)
  write(150,rec=ik,err=10) vkl(:,ik),nstsv,pmat
  close(150)
  exit
10 continue
  if (i.eq.2) then
    write(*,*)
    write(*,'("Error(putpmat): unable to write to PMAT.OUT")')
    write(*,*)
    stop
  end if
  close(150)
end do
deallocate(pmat)
return
end subroutine

