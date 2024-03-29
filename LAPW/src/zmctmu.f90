SUBROUTINE zmctmu(tcr,l,n,a,b,ld,c)
  IMPLICIT NONE 
  ! arguments
  logical, intent(in) :: tcr
  INTEGER, intent(in) :: l,n
  COMPLEX(8), intent(in) :: a(l,n),b(l,n)
  INTEGER, intent(in) :: ld
  COMPLEX(8), intent(out) :: c(*)
  ! local variables
  INTEGER :: l2,i,j,k
  ! external functions
  REAL(8) ddot
  COMPLEX(8) zdotc
  external ddot,zdotc

  IF(tcr) THEN 
    ! matrix c is real valued
    l2=2*l
    DO j=1,n
      k=(j-1)*ld
      DO i=1,j
        k=k+1
        c(k)=c(k)+ddot(l2,a(:,i),1,b(:,j),1)
      ENDDO 
    ENDDO 
  ELSE 
    ! matrix c is complex valued
    DO j=1,n
      k=(j-1)*ld
      DO i=1,j
        k=k+1
        c(k)=c(k)+zdotc(l,a(:,i),1,b(:,j),1)
      ENDDO 
    ENDDO 
  ENDIF 
  RETURN 
END SUBROUTINE 

