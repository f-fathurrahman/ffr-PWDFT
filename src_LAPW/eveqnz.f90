SUBROUTINE eveqnz(n,ld,a,w)
  use modmain
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: n,ld
  COMPLEX(8), intent(inout) :: a(ld,n)
  REAL(8), intent(out) :: w(n)
  ! local variables
  INTEGER :: liwork,lrwork
  INTEGER :: lwork,info,nthd
  ! ALLOCATABLE arrays
  INTEGER, ALLOCATABLE :: iwork(:)
  REAL(8), ALLOCATABLE :: rwork(:)
  COMPLEX(8), ALLOCATABLE :: work(:)
  
  SELECT CASE(evtype)
  CASE(0)
    ! use the LAPACK routine zheev
    lwork=2*n
    ALLOCATE(rwork(3*n),work(lwork))
    ! enable MKL parallelism
    CALL zheev('V','U',n,a,ld,w,work,lwork,rwork,info)
    IF(info.ne.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(eveqnz): diagonalisation failed")')
      WRITE(*,'(" ZHEEV RETURN ed INFO = ",I8)') info
      WRITE(*,*)
      STOP 
    ENDIF 
    DEALLOCATE(rwork,work)
  CASE(1)
    ! use the divide-and-conquer LAPACK routine zheevd
    liwork=5*n+3
    lrwork=2*n**2+5*n+1
    lwork=n**2+2*n
    ALLOCATE(iwork(liwork),rwork(lrwork),work(lwork))
    ! enable MKL parallelism
    CALL zheevd('V','U',n,a,ld,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    IF(info.ne.0) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(eveqnz): diagonalisation failed")')
      WRITE(*,'(" ZHEEVD RETURN ed INFO = ",I8)') info
      WRITE(*,*)
      stop
    ENDIF 
    DEALLOCATE(iwork,rwork,work)
  CASE DEFAULT 
    WRITE(*,*)
    WRITE(*,'("Error(eveqnz): evtype not defined : ",I8)') evtype
    WRITE(*,*)
    STOP 
  END SELECT 
  RETURN 
END SUBROUTINE 

