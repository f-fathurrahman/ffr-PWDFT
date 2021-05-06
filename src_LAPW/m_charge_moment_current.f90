module m_charge_moment_current

use m_atoms, only: maxspecies

!----------------------------------------------!
!     charge, moment and current variables     !
!----------------------------------------------!
! tolerance for error in total charge
real(8) epschg
! total nuclear charge
real(8) chgzn
! core charges
real(8) chgcr(maxspecies)
! total core charge
real(8) chgcrtot
! core leakage charge
real(8), allocatable :: chgcrlk(:)
! total valence charge
real(8) chgval
! excess charge
real(8) chgexs
! total charge
real(8) chgtot
! calculated total charge
real(8) chgcalc
! interstitial region charge
real(8) chgir
! muffin-tin charges
real(8), allocatable :: chgmt(:)
! total muffin-tin charge
real(8) chgmttot
! effective Wigner radius
real(8) rwigner
! total moment
real(8) momtot(3)
! total moment magnitude
real(8) momtotm
! interstitial region moment
real(8) momir(3)
! muffin-tin moments
real(8), allocatable :: mommt(:,:)
! total muffin-tin moment
real(8) mommttot(3)
! total current
real(8) curtot(3)
! total current magnitude
real(8) curtotm

end module 
