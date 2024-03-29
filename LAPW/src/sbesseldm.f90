! !INPUT/OUTPUT PARAMETERS:
!   m    : order of derivatve (in,integer)
!   lmax : maximum order of Bessel function (in,integer)
!   x    : real argument (in,real)
!   djl  : array of RETURN ed values (out,real(0:lmax))
SUBROUTINE sbesseldm(m,lmax,x,djl)
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: m,lmax
  REAL(8), intent(in) :: x
  REAL(8), intent(out) :: djl(0:lmax)
  ! local variables
  INTEGER, parameter :: maxm=6,maxns=20
  INTEGER i,j,l,i0
  REAL(8) t1,sum,x2
  INTEGER a(0:maxm+1),a1(0:maxm+1)
  INTEGER b(0:maxm+1),b1(0:maxm+1)
  ! automatic arrays
  REAL(8) jl(0:lmax+1)
  ! external functions
  REAL(8) factnm,factr
  external factnm,factr

  IF( (m < 0) .or. (m > maxm)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(sbesseldm): m out of range : ",I8)') m
    WRITE(*,*)
    STOP 
  ENDIF 

  IF((lmax < 0) .or. (lmax > 30)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(sbesseldm): lmax out of range : ",I8)') lmax
    WRITE(*,*)
    STOP
  ENDIF 
  
  IF((x < 0.d0) .or. (x > 1.d5)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(sbesseldm): x out of range : ",G18.10)') x
    WRITE(*,*)
    STOP 
  ENDIF 
  
  IF(m==0) THEN 
    CALL sbessel(lmax,x,djl)
    RETURN 
  ENDIF 

  IF(x > 1.d0) THEN 
    CALL sbessel(lmax+1,x,jl)
    DO l=0,lmax
      a(1:m+1)=0
      a(0)=1
      a1(0:m+1)=0
      DO i=1,m
        b(0)=0
        b1(0)=0
        DO j=0,i
          b(j+1)=a(j)*(l-j)
          b1(j+1)=-a1(j)*(j+l+2)
        ENDDO 
        DO j=0,i
          b1(j)=b1(j)-a(j)
          b(j)=b(j)+a1(j)
        ENDDO 
        a(0:i+1)=b(0:i+1)
        a1(0:i+1)=b1(0:i+1)
      ENDDO 
      t1=1.d0
      sum=dble(a(0))*jl(l)+dble(a1(0))*jl(l+1)
      DO i=1,m+1
        t1=t1*x
        sum=sum+(dble(a(i))*jl(l)+dble(a1(i))*jl(l+1))/t1
      ENDDO 
      djl(l)=sum
    ENDDO 
  ELSE 
    x2=x**2
    DO l=0,lmax
      i0=max((m-l+1)/2,0)
      j=2*i0+l-m
      IF(j == 0) THEN 
        t1=1.d0
      else
        t1=x**j
      ENDIF 
      t1=factr(j+m,j)*t1/(factnm(i0,1)*factnm(j+l+m+1,2)*dble((-2)**i0))
      sum=t1
      DO i=i0+1,maxns
        j=2*i+l
        t1=-t1*dble((j-1)*j)*x2/dble((j-l)*(j-m-1)*(j-m)*(j+l+1))
        IF(abs(t1) <= 1.d-40) GOTO 10
        sum=sum+t1
      ENDDO 
      10 CONTINUE 
      djl(l)=sum
    ENDDO 
  ENDIF 
  RETURN 
END SUBROUTINE 
