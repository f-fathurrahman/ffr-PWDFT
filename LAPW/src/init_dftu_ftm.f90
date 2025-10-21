subroutine init_dftu_ftm()

use modmain
use moddftu
implicit none

integer :: i, ia, ias, is, l

!-------------------------------------------------!
!     DFT+U and fixed tensor moment variables     !
!-------------------------------------------------!
if ((dftu.ne.0).or.(ftmtype.ne.0)) then
! density matrix elements in each muffin-tin
  if (allocated(dmatmt)) deallocate(dmatmt)
  allocate(dmatmt(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
! potential matrix elements in each muffin-tin
  if (allocated(vmatmt)) deallocate(vmatmt)
  allocate(vmatmt(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
! zero the potential matrix
  vmatmt(:,:,:,:,:)=0.d0
! require the potential matrix elements be calculated
  tvmatmt=.true.
! flags for non-zero muffin-tin potential matrices
  if (allocated(tvmmt)) deallocate(tvmmt)
  allocate(tvmmt(0:lmaxdm,natmtot))
  tvmmt(:,:)=.false.
! require second-variational eigenvectors
  tevecsv=.true.
end if
if (dftu.ne.0) then
! DFT+U energy for each atom
  if (allocated(engyadu)) deallocate(engyadu)
  allocate(engyadu(natmmax,ndftu))
! interpolation constants (alpha)
  if (allocated(alphadu)) deallocate(alphadu)
  allocate(alphadu(natmmax,ndftu))
! flag the muffin-tin potential matrices which are non-zero
  do i=1,ndftu
    is=idftu(1,i)
    if (is.gt.nspecies) then
      write(*,*)
      write(*,'("Error(init0): invalid species number : ",I8)') is
      write(*,*)
      stop
    end if
    l=idftu(2,i)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      tvmmt(l,ias)=.true.
    end do
  end do
end if
if (ftmtype.ne.0) then
! allocate and zero the fixed tensor moment potential array
  if (allocated(vmftm)) deallocate(vmftm)
  allocate(vmftm(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
  vmftm(:,:,:,:,:)=0.d0
! flag the muffin-tin potential matrices which are non-zero
  do i=1,ntmfix
    is=itmfix(1,i)
    ia=itmfix(2,i)
    ias=idxas(ia,is)
    l=itmfix(3,i)
    tvmmt(l,ias)=.true.
  end do
end if

end subroutine