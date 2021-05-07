! !DESCRIPTION:
!   Interface to the grid requirements of the fast Fourier transform routine.
!   Most routines restrict $n$ to specific prime factorisations. This routine
!   RETURN s the next largest grid size allowed by the FFT routine.
SUBROUTINE nfftifc(n)
! !INPUT/OUTPUT PARAMETERS:
!   n : required/avalable grid size (inout,integer)
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(inout) :: n
  ! local variables
  INTEGER :: i,j
  ! currently we use primes 2, 3 and 5
  INTEGER, parameter :: np=3
  INTEGER p(np)
  data p / 2,3,5 /

  IF(n <= 0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(nfftifc): n <= 0 : ",I8)') n
    WRITE(*,*)
    STOP 
  ENDIF 
  10 CONTINUE 
  
  i=n
  DO j=1,np
    DO WHILE( mod(i,p(j))==0 )
      i = i/p(j)
    ENDDO 
  ENDDO 
  
  IF(i /= 1) THEN 
    n = n + 1
    GOTO 10
  ENDIF 

END SUBROUTINE 
