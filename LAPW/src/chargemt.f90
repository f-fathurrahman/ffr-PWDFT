
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine chargemt
use modmain
implicit none
! local variables
integer is,ias,nr,nri
real(8) t1
! automatic arrays
real(8) fr(nrmtmax)
chgmttot=0.d0
do ias=1,natmtot
  is = idxis(ias)
  nr = nrmt(is)
  nri = nrmti(is)
  ! extract the l=m=0 component from the muffin-tin density
  call rfmtlm(1, nr, nri, rhomt(:,ias), fr)
  ! integrate to the muffin-tin radius
  t1 = dot_product(wrmt(1:nr,is), fr(1:nr))
  chgmt(ias) = fourpi*y00*t1
  write(*,'(1x,A,I4,A,F18.10)') 'ias = ', ias, 'chgmt(ias) = ', chgmt(ias)
  chgmttot = chgmttot + chgmt(ias)
enddo
return
end subroutine

