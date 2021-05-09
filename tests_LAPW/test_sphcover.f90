PROGRAM main
  IMPLICIT NONE 
  INTEGER, PARAMETER :: n=100
  INTEGER :: i
  REAL(8) :: tp(2,n)

  CALL sphcover(n, tp)
  DO i=1,n
    WRITE(*,'(1x,I8,2F18.10)') i, tp(1,i), tp(2,i)
  ENDDO 

END PROGRAM 
