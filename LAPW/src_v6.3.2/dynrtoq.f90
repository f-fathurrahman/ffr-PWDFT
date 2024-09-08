
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dynrtoq(vpl,dynr,dynp)
use modmain
use modphonon
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(in) :: dynr(nbph,nbph,nqptnr)
complex(8), intent(out) :: dynp(nbph,nbph)
! local variables
integer i1,i2,i3,ir,i,j
real(8) t1
complex(8) z1
dynp(:,:)=0.d0
! loop over R-vectors
ir=0
do i3=ngridq(3)/2-ngridq(3)+1,ngridq(3)/2
  do i2=ngridq(2)/2-ngridq(2)+1,ngridq(2)/2
    do i1=ngridq(1)/2-ngridq(1)+1,ngridq(1)/2
      ir=ir+1
      t1=-twopi*(vpl(1)*dble(i1)+vpl(2)*dble(i2)+vpl(3)*dble(i3))
      z1=cmplx(cos(t1),sin(t1),8)
      do i=1,nbph
        do j=1,nbph
          dynp(i,j)=dynp(i,j)+z1*dynr(i,j,ir)
        end do
      end do
    end do
  end do
end do
! symmetrise the dynamical matrix
call dynsym(vpl,dynp)
return
end subroutine

