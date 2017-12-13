! adapted from: http://fortranwiki.org/fortran/show/Matrix+inversion

SUBROUTINE inv_m3x3(A, B)
  IMPLICIT NONE 
  REAL(8) :: A(3,3)   !! Matrix
  REAL(8) :: B(3,3)   !! Inverse matrix
  REAL(8) :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1.d0/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
               - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
               + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

END SUBROUTINE 

