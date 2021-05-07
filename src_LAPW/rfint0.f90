
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfint0(rf0,rfmt,rfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rf0
real(8), intent(inout) :: rfmt(npmtmax,natmtot),rfir(ngtot)
! local variables
integer is,ias
integer nr,nri,ir,i
real(8) t1
! external functions
real(8) rfint
external rfint
t1=rfint(rfmt,rfir)
t1=rf0-t1/omega
rfir(:)=rfir(:)+t1
t1=t1/y00
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  i=1
  do ir=1,nri
    rfmt(i,ias)=rfmt(i,ias)+t1
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    rfmt(i,ias)=rfmt(i,ias)+t1
    i=i+lmmaxo
  end do
end do
return
end subroutine

