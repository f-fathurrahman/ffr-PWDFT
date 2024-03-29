! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   rflm : coefficients of real spherical harmonic expansion
!          (in,real((lmax+1)**2)))
!   zflm : coefficients of complex spherical harmonic expansion
!          (out,complex((lmax+1)**2)))
PURE SUBROUTINE rtozflm(lmax,rflm,zflm)
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: lmax
  REAL(8), intent(in) :: rflm(*)
  COMPLEX(8), intent(out) :: zflm(*)
  ! local variables
  INTEGER l,m,lm1,lm2
  ! real constant 1/sqrt(2)
  REAL(8), parameter :: c1=0.7071067811865475244d0
  lm1=0
  DO l=0,lmax
    lm2=lm1+2*(l+1)
    DO m=-l,-1
      lm1=lm1+1
      lm2=lm2-1
      IF(mod(m,2).ne.0) THEN 
        zflm(lm1)=c1*cmplx(-rflm(lm2),-rflm(lm1),8)
      else
        zflm(lm1)=c1*cmplx(rflm(lm2),-rflm(lm1),8)
      ENDIF 
    ENDDO 
    lm1=lm1+1
    lm2=lm2-1
    zflm(lm1)=rflm(lm1)
    DO m=1,l
      lm1=lm1+1
      lm2=lm2-1
      IF(mod(m,2).ne.0) THEN 
        zflm(lm1)=c1*cmplx(rflm(lm1),-rflm(lm2),8)
      else
        zflm(lm1)=c1*cmplx(rflm(lm1),rflm(lm2),8)
      ENDIF 
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 