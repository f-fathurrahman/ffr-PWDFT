
! Copyright (C) 2003-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rfmtinp
! !INTERFACE:
real(8) function rfmtinp(nr,nri,wr,rfmt1,rfmt2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr    : number of radial mesh points (in,integer)
!   nri   : number of radial mesh points on the inner part of the muffin-tin
!           (in,integer)
!   wr    : weights for integration on radial mesh (in,real(nr))
!   rfmt1 : first real function inside muffin-tin (in,real(*))
!   rfmt2 : second real function inside muffin-tin (in,real(*))
! !DESCRIPTION:
!   Calculates the inner product of two real functions in the muffin-tin. So
!   given two real functions of the form
!   $$ f({\bf r})=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}f_{lm}(r)R_{lm}
!    (\hat{\bf r}) $$
!   where $R_{lm}$ are the real spherical harmonics, the function returns
!   $$ I=\int\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}f_{lm}^1(r)f_{lm}^2(r)r^2
!    dr\;. $$
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr)
real(8), intent(in) :: rfmt1(*),rfmt2(*)
! local variables
integer n,ir,i
! automatic arrays
real(8) fr(nr)
i=1
! inner part of muffin-tin
if (lmaxi.eq.1) then
  do ir=1,nri
    fr(ir)=rfmt1(i)*rfmt2(i) &
          +rfmt1(i+1)*rfmt2(i+1) &
          +rfmt1(i+2)*rfmt2(i+2) &
          +rfmt1(i+3)*rfmt2(i+3)
    i=i+4
  end do
else
  n=lmmaxi-1
  do ir=1,nri
    fr(ir)=dot_product(rfmt1(i:i+n),rfmt2(i:i+n))
    i=i+lmmaxi
  end do
end if
! outer part of muffin-tin
n=lmmaxo-1
do ir=nri+1,nr
  fr(ir)=dot_product(rfmt1(i:i+n),rfmt2(i:i+n))
  i=i+lmmaxo
end do
! integrate
rfmtinp=dot_product(wr(:),fr(:))
return
end function
!EOC

