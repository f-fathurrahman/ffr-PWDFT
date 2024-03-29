! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   zflm : coefficients of complex spherical harmonic expansion
!          (in,complex((lmax+1)**2)))
!   rflm : coefficients of real spherical harmonic expansion
!          (out,real((lmax+1)**2)))
PURE SUBROUTINE z_to_rf_lm(lmax, zflm, rflm)
  IMPLICIT NONE 
  ! arguments
  INTEGER, INTENT(in) :: lmax
  COMPLEX(8), INTENT(in) :: zflm(*)
  REAL(8), INTENT(out) :: rflm(*)
  ! local variables
  INTEGER :: l,m,lm1,lm2
  ! real constant 1/sqrt(2)
  REAL(8), PARAMETER :: c1=0.7071067811865475244d0
  
  lm1 = 0
  DO l=0,lmax
    lm2 = lm1 + 2*(l+1)
    DO m=-l,-1
      lm1 = lm1+1
      lm2 = lm2-1
      IF(mod(m,2) /= 0) THEN 
        rflm(lm1) = -c1*(aimag(zflm(lm1)) + aimag(zflm(lm2)))
      ELSE 
        rflm(lm1) = c1*(aimag(zflm(lm2)) - aimag(zflm(lm1)))
      ENDIF 
    ENDDO 
    lm1 = lm1+1
    lm2 = lm2-1
    rflm(lm1) = dble(zflm(lm1))
    DO m = 1,l
      lm1=lm1+1
      lm2=lm2-1
      IF( mod(m,2) /= 0) THEN 
        rflm(lm1) = c1*(dble(zflm(lm1)) - dble(zflm(lm2)))
      ELSE 
        rflm(lm1) = c1*(dble(zflm(lm1)) + dble(zflm(lm2)))
      ENDIF 
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 