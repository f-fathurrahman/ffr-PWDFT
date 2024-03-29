! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   v    : input vector (in,real(3))
!   rlm  : array of real spherical harmonics (out,real((lmax+1)**2))
SUBROUTINE genrlmv(lmax,v,rlm)
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: lmax
  REAL(8), intent(in) :: v(3)
  REAL(8), intent(out) :: rlm(*)
  ! local variables
  INTEGER l,m,lm
  REAL(8), parameter :: sqtwo=1.4142135623730950488d0
  ! automatic arrays
  COMPLEX(8) ylm((lmax+1)**2)
  IF((lmax.lt.0).or.(lmax.gt.50)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(genrlmv): lmax out of range : ",I8)') lmax
    WRITE(*,*)
    STOP 
  ENDIF 
  ! generate complex spherical harmonics
  CALL genylmv(lmax,v,ylm)
  ! convert to real spherical harmonics
  lm=0
  DO l=0,lmax
    DO m=-l,-1
      lm=lm+1
      rlm(lm)=sqtwo*aimag(ylm(lm))
    ENDDO 
    lm=lm+1
    rlm(lm)=dble(ylm(lm))
    DO m=1,l
      lm=lm+1
      rlm(lm)=sqtwo*dble(ylm(lm))
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 
