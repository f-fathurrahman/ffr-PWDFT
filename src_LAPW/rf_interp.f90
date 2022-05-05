! !INPUT/OUTPUT PARAMETERS:
!   ni : number of input points (in,integer)
!   xi : input abscissa array (in,real(ni))
!   fi : input data array (in,real(ni)
!   no : number of output points (in,integer)
!   xo : output abscissa array (in,real(ni))
!   fo : output interpolated function (out,real(no))
SUBROUTINE rf_interp(ni,xi,fi,no,xo,fo)
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ni
  REAL(8), intent(in) :: xi(ni),fi(ni)
  INTEGER, intent(in) :: no
  REAL(8), intent(in) :: xo(no)
  REAL(8), intent(out) :: fo(no)
  ! local variables
  INTEGER :: i,j,k,l
  REAL(8) :: x,dx
  ! automatic arrays
  REAL(8) :: cf(3,ni)
  IF(ni <= 0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(rfinterp): invalid number of input points : ",I8)') ni
    WRITE(*,*)
    STOP 
  ENDIF 
  IF(no <= 0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(rfinterp): invalid number of output points : ",I8)') no
    WRITE(*,*)
    STOP 
  ENDIF 
  IF(ni==1) THEN 
    fo(:)=fi(1)
    RETURN 
  ENDIF 
  ! compute the spline coefficients
  CALL spline(ni,xi,fi,cf)
  ! evaluate spline at output points
  i=1
  DO l=1,no
    x=xo(l)
    IF(i >= ni) i=1
    IF(x < xi(i)) GOTO 10
    IF(x <= xi(i+1)) GOTO 30
    ! binary search
    10 CONTINUE 
    i=1
    j=ni
    20 CONTINUE 
    k=(i+j)/2
    IF(x < xi(k)) THEN 
      j=k
    ELSE 
      i=k
    ENDIF 
    IF(j > i+1) goto 20
    30 CONTINUE 
    dx=x-xi(i)
    fo(l)=fi(i)+dx*(cf(1,i)+dx*(cf(2,i)+dx*cf(3,i)))
  ENDDO 
  RETURN 
END SUBROUTINE 
