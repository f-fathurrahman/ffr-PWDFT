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
  ! find the interstitial charge
  t1=dot_product(rhoir(:),cfunir(:))
  chgir=t1*omega/dble(ngtot)
  ! total calculated charge
  chgcalc=chgmttot+chgir

  WRITE(*,*) 'calculated total charge', chgcalc
  RETURN 
END SUBROUTINE 
