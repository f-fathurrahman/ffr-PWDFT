PROGRAM test_mod
  IMPLICIT NONE 
  REAL(8) :: pi, p, dp

  pi = 4.d0*atan(1.d0)
  dp = pi*(1.0 - sqrt(5.d0))
  WRITE(*,*) 'dp = ', dp

  p = 0.d0
  p = p + dp
  WRITE(*,*) 'p = ', p
  WRITE(*,*) 'mod(p,2pi) = ', mod(p, 2.d0*pi)
END PROGRAM 

