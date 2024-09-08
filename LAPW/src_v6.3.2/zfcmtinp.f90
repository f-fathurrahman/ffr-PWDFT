
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function zfcmtinp(nr,nri,wr,zfmt1,zfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr)
complex(8), intent(in) :: zfmt1(*),zfmt2(*)
! local variables
integer ir,i
real(8) t1
complex(8) z1
! automatic arrays
real(8) fr1(nr),fr2(nr)
! external functions
complex(8) zdotc
external zdotc
! compute the dot-products for each radial point
i=1
if (lmaxi.eq.1) then
  do ir=1,nri
    z1=conjg(zfmt1(i))*zfmt2(i) &
      +conjg(zfmt1(i+1))*zfmt2(i+1) &
      +conjg(zfmt1(i+2))*zfmt2(i+2) &
      +conjg(zfmt1(i+3))*zfmt2(i+3)
    fr1(ir)=dble(z1); fr2(ir)=aimag(z1)
    i=i+4
  end do
else
  do ir=1,nri
    z1=zdotc(lmmaxi,zfmt1(i),1,zfmt2(i),1)
    fr1(ir)=dble(z1); fr2(ir)=aimag(z1)
    i=i+lmmaxi
  end do
end if
t1=dble(lmmaxi)/dble(lmmaxo)
do ir=nri+1,nr
  z1=t1*zdotc(lmmaxo,zfmt1(i),1,zfmt2(i),1)
  fr1(ir)=dble(z1); fr2(ir)=aimag(z1)
  i=i+lmmaxo
end do
! integrate over r
t1=fourpi/dble(lmmaxi)
zfcmtinp=t1*cmplx(dot_product(wr(:),fr1(:)),dot_product(wr(:),fr2(:)),8)
return
end function

