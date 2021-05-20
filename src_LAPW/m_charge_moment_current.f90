MODULE m_charge_moment_current

USE m_atoms, only: maxspecies

!----------------------------------------------!
!     charge, moment and current variables     !
!----------------------------------------------!
! tolerance for error in total charge
REAL(8) epschg
! total nuclear charge
REAL(8) chgzn
! core charges
REAL(8) chgcr(maxspecies)
! total core charge
REAL(8) chgcrtot
! core leakage charge
REAL(8), ALLOCATABLE :: chgcrlk(:)
! total valence charge
REAL(8) chgval
! excess charge
REAL(8) chgexs
! total charge
REAL(8) chgtot
! calculated total charge
REAL(8) chgcalc
! interstitial region charge
REAL(8) chgir
! muffin-tin charges
REAL(8), ALLOCATABLE :: chgmt(:)
! total muffin-tin charge
REAL(8) chgmttot
! effective Wigner radius
REAL(8) rwigner
! total moment
REAL(8) momtot(3)
! total moment magnitude
REAL(8) momtotm
! interstitial region moment
REAL(8) momir(3)
! muffin-tin moments
REAL(8), ALLOCATABLE :: mommt(:,:)
! total muffin-tin moment
REAL(8) mommttot(3)
! total current
REAL(8) curtot(3)
! total current magnitude
REAL(8) curtotm

end module 
