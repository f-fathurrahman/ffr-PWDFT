
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmat(tspndg,tlmdg,lmin,lmax,ld,dmat)
use modmain
implicit none
! arguments
logical, intent(in) :: tspndg,tlmdg
integer, intent(in) :: lmin,lmax,ld
complex(8), intent(out) :: dmat(ld,nspinor,ld,nspinor,natmtot)
! local variables
integer ik,ispn,ist
integer ias
real(8) wo
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:),evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: dmatk(:,:,:,:,:)

! zero the density matrix
dmat(:,:,:,:,:)=0.d0

! begin parallel loop over k-points

allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
allocate(dmatk(ld,nspinor,ld,nspinor,nstsv))

do ik=1,nkpt

! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),vgkc(:,:,ispn,ik),gkc(:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! get the eigenvectors from file
  call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
! begin loop over atoms and species
  do ias=1,natmtot
    call gendmatk(tspndg,tlmdg,lmin,lmax,ias,ngk(:,ik),apwalm,evecfv,evecsv, &
     ld,dmatk)
    do ist=1,nstsv
      wo=wkpt(ik)*occsv(ist,ik)
      if (wo.lt.epsocc) cycle
      dmat(:,:,:,:,ias)=dmat(:,:,:,:,ias)+wo*dmatk(:,:,:,:,ist)
    end do
  end do
end do
deallocate(apwalm,evecfv,evecsv,dmatk)

!ffr: mpi stuffs removed

! symmetrise the density matrix
call symdmat(lmax,ld,dmat)

return
end subroutine

