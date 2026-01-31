SUBROUTINE charge()
  ! !USES:
  use modmain
  ! !DESCRIPTION:
  !   Computes the muffin-tin, interstitial and total charges by integrating the
  !   density.
  IMPLICIT NONE 
  ! local variables
  REAL(8) :: t1
  
  ! find the muffin-tin charges
  CALL chargemt()
  write(*,*) 'chgmttot = ', chgmttot
  
  ! find the interstitial charge
  write(*,*) 'sum(rhoir) = ', sum(rhoir)
  t1 = dot_product(rhoir(:),cfunir(:))
  write(*,*) 'for interstitial charge: t1 = ', t1
  chgir = t1*omega/dble(ngtot)
  write(*,*) 'chgir = ', chgir

  ! total calculated charge
  chgcalc = chgmttot + chgir

  WRITE(*,*) 'calculated total charge', chgcalc
  RETURN 
END SUBROUTINE 
