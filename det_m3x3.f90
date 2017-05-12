! adapted from: http://fortranwiki.org/fortran/show/Matrix+inversion

FUNCTION det_m3x3(A) RESULT(det)
  IMPLICIT NONE 
  REAL(8) :: A(3,3)   !! Matrix
  REAL(8) :: det

  det = (  A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
         - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
         + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1) )

END FUNCTION 

