
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmdedc
! !INTERFACE:
subroutine rdmdedc(dedc)
! !USES:
use modmain
use modrdm
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   dedc : energy derivative (out,complex(nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Calculates the derivative of the total energy w.r.t. the second-variational
!   coefficients {\tt evecsv}.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(out) :: dedc(nstsv,nstsv,nkpt)
! local variables
integer ik,ist,nthd
! allocatable arrays
complex(8), allocatable :: evecsv(:,:),c(:,:)
call holdthd(nkpt,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecsv,c,ist) &
!$OMP NUM_THREADS(nthd)
allocate(evecsv(nstsv,nstsv),c(nstsv,nstsv))
!$OMP DO
do ik=1,nkpt
! get the eigenvectors from file
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
! kinetic and Coulomb potential contribution
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,vclmat(:,:,ik),nstsv, &
   zzero,c,nstsv)
  do ist=1,nstsv
    dedc(:,ist,ik)=occsv(ist,ik)*(dkdc(:,ist,ik)+c(:,ist))
  end do
end do
!$OMP END DO
deallocate(evecsv,c)
!$OMP END PARALLEL
call freethd(nthd)
! exchange-correlation contribution
call rdmdexcdc(dedc)
return
end subroutine
!EOC

