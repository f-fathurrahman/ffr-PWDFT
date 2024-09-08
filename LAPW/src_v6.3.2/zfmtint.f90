
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function zfmtint(nr,nri,wr,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr)
complex(8), intent(in) :: zfmt(*)
! local variables
integer ir,i
complex(8) z1
! automatic arrays
real(8) fr1(nr),fr2(nr)
i=1
do ir=1,nri
  z1=zfmt(i)
  fr1(ir)=dble(z1); fr2(ir)=aimag(z1)
  i=i+lmmaxi
end do
do ir=nri+1,nr
  z1=zfmt(i)
  fr1(ir)=dble(z1); fr2(ir)=aimag(z1)
  i=i+lmmaxo
end do
! integrate over r
z1=cmplx(dot_product(wr(:),fr1(:)),dot_product(wr(:),fr2(:)),8)
zfmtint=fourpi*y00*z1
return
end function

