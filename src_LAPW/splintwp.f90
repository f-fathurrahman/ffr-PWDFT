PURE SUBROUTINE splintwp(n,wp,f,g)
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: n
  REAL(8), intent(in) :: wp(4,n),f(n)
  REAL(8), intent(out) :: g(n)
  ! local variables
  INTEGER i
  g(1)=0.d0
  g(2)=wp(1,2)*f(1)+wp(2,2)*f(2)+wp(3,2)*f(3)+wp(4,2)*f(4)
  DO i=3,n-1
    g(i)=g(i-1)+wp(1,i)*f(i-2)+wp(2,i)*f(i-1)+wp(3,i)*f(i)+wp(4,i)*f(i+1)
  ENDDO 
  g(n)=g(n-1)+wp(1,n)*f(n-3)+wp(2,n)*f(n-2)+wp(3,n)*f(n-1)+wp(4,n)*f(n)
  RETURN 
END SUBROUTINE 

