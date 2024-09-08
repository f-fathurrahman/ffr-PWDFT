
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine orthevsv(evecsv)
use modmain
implicit none
! arguments
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
integer ist,jst
real(8) t1
complex(8) z1
! external functions
complex(8) zdotc
external zdotc
! perform Gram-Schmidt orthonormalisation
do ist=1,nstsv
  do jst=1,ist-1
    z1=-zdotc(nstsv,evecsv(:,jst),1,evecsv(:,ist),1)
    call zaxpy(nstsv,z1,evecsv(:,jst),1,evecsv(:,ist),1)
  end do
  t1=dble(zdotc(nstsv,evecsv(:,ist),1,evecsv(:,ist),1))
  if (t1.gt.0.d0) then
    t1=1.d0/sqrt(t1)
    call zdscal(nstsv,t1,evecsv(:,ist),1)
  end if
end do
return
end subroutine

