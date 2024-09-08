
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function rfinpc(ld,rfmt1,rfir1,rfmt2,rfir2)
use modmain
use modomp
implicit none
integer, intent(in) :: ld
real(8), intent(in) :: rfmt1(ld,natmtot),rfir1(ngtot)
real(8), intent(in) :: rfmt2(ld,natmtot),rfir2(ngtot)
! local variables
integer is,ias,ir,nthd
! external functions
real(8) rfmtinp
external rfmtinp
! interstitial contribution
rfinpc=0.d0
do ir=1,ngtot
  rfinpc=rfinpc+rfir1(ir)*rfir2(ir)*cfunir(ir)
end do
rfinpc=rfinpc*omega/dble(ngtot)
! muffin-tin contribution
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is) REDUCTION(+:rfinpc) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  rfinpc=rfinpc+rfmtinp(nrcmt(is),nrcmti(is),wrcmt(:,is),rfmt1(:,ias), &
   rfmt2(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
return
end function

