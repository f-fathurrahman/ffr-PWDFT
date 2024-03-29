! !DESCRIPTION:
!   Finds the permutation index {\tt idx} which sorts the real array {\tt x}
!   into ascending order. No sorting of the array {\tt x} itself is performed.
!   Uses the heapsort algorthim.
SUBROUTINE sortidx(n,x,idx)
! !INPUT/OUTPUT PARAMETERS:
!   n   : number of elements in array (in,integer)
!   x   : real array (in,real(n))
!   idx : permutation index (out,integer(n))
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: n
  REAL(8), intent(in) :: x(n)
  INTEGER, intent(out) :: idx(n)
  ! local variables
  INTEGER i,j,k,l,m
  ! tolerance for deciding if one number is smaller than another
  REAL(8), parameter :: eps=1.d-14
  
  IF(n <= 0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(sortidx): n <= 0 : ",I8)') n
    WRITE(*,*)
    STOP 
  ENDIF 
  DO i=1,n
    idx(i)=i
  ENDDO 
  IF(n==1) RETURN 
  l = n/2 + 1
  k = n
  
  10 continue
  
  IF(l > 1) THEN 
    l=l-1
    m=idx(l)
  ELSE 
    m=idx(k)
    idx(k)=idx(1)
    k=k-1
    IF(k == 1) THEN 
      idx(1) = m
      RETURN 
    ENDIF 
  ENDIF 

  i = l
  j = l + l
  
  20 continue
  
  IF(j <= k) THEN 
    IF(j < k) THEN 
      IF(x(idx(j)) < x(idx(j+1))+eps) j=j+1
    ENDIF 
    IF(x(m) < x(idx(j))+eps) THEN 
      idx(i) = idx(j)
      i=j
      j=j+j
    ELSE 
      j=k+1
    ENDIF 
    GOTO 20
  ENDIF 
  idx(i) = m
  GOTO 10
END SUBROUTINE 

