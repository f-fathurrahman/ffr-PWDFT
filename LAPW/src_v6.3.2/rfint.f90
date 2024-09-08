
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function rfint(rfmt,rfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
! local variables
integer is,ias,nr,nri
real(8) t1
! automatic arrays
real(8) fr(nrmtmax)
! interstitial contribution
rfint=dot_product(rfir(:),cfunir(:))
rfint=rfint*omega/dble(ngtot)
! muffin-tin contribution
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
! extract the l=m=0 component
  call rfmtlm(1,nr,nri,rfmt(:,ias),fr)
! integrate to the muffin-tin radius
  t1=dot_product(wrmt(1:nr,is),fr(1:nr))
  rfint=rfint+fourpi*y00*t1
end do
return
end function

