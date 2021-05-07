! !INPUT/OUTPUT PARAMETERS:
!   n  : number of required points (in,integer)
!   tp : (theta, phi) coordinates (out,real(2,n))
SUBROUTINE sphcover(n,tp)
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: n
  REAL(8), intent(out) :: tp(2,n)
  ! local variables
  INTEGER k
  REAL(8), parameter :: pi=3.1415926535897932385d0
  REAL(8) z,dz,p,dp
  IF(n.le.0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(sphcover): n <= 0 : ",I8)') n
    WRITE(*,*)
    stop
  ENDIF 
  dz=2.d0/dble(n)
  z=1.d0-dz/2.d0
  tp(1,1)=acos(z)
  dp=pi*(1.d0-sqrt(5.d0))
  p=0.d0
  tp(2,1)=p
  DO k=2,n
    z=z-dz
    tp(1,k)=acos(z)
    p=p+dp
    tp(2,k)=mod(p,2.d0*pi)
  ENDDO 
  RETURN 
END SUBROUTINE 