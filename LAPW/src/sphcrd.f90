! !INPUT/OUTPUT PARAMETERS:
!   v  : input vector (in,real(3))
!   r  : length of v (out,real)
!   tp : (theta, phi) coordinates (out,real(2))
PURE SUBROUTINE sphcrd(v,r,tp)
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(in) :: v(3)
  REAL(8), intent(out) :: r,tp(2)
  ! local variables
  REAL(8), parameter :: eps=1.d-14
  REAL(8) t1
  r=sqrt(v(1)**2+v(2)**2+v(3)**2)
  IF(r.gt.eps) THEN 
    t1=v(3)/r
    IF(t1.ge.1.d0) THEN 
      tp(1)=0.d0
    ELSEIF(t1.le.-1.d0) THEN 
      tp(1)=3.1415926535897932385d0
    else
      tp(1)=acos(t1)
    ENDIF 
    IF((abs(v(1)).gt.eps).or.(abs(v(2)).gt.eps)) THEN 
      tp(2)=atan2(v(2),v(1))
    else
      tp(2)=0.d0
    ENDIF 
  else
    tp(1)=0.d0
    tp(2)=0.d0
  ENDIF 
  RETURN 
END SUBROUTINE 