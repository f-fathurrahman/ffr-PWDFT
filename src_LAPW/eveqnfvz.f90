SUBROUTINE eveqnfvz(nmatp,h,o,evalfv,evecfv)
  USE m_states, ONLY: nstfv, evaltol
  USE m_timing, ONLY: timefv
  USE m_hamiltonian, ONLY: nmatmax
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nmatp
  COMPLEX(8), intent(in) :: h(*),o(*)
  REAL(8), intent(out) :: evalfv(nstfv)
  COMPLEX(8), intent(out) :: evecfv(nmatmax,nstfv)
  ! local variables
  INTEGER :: i,m
  INTEGER :: lwork,info
  REAL(8) :: vl,vu
  REAL(8) :: ts0,ts1
  ! ALLOCATABLE arrays
  INTEGER, ALLOCATABLE :: iwork(:),ifail(:)
  REAL(8), ALLOCATABLE :: w(:),rwork(:)
  COMPLEX(8), ALLOCATABLE :: work(:)

  CALL timesec(ts0)

  ALLOCATE(iwork(5*nmatp),ifail(nmatp))
  ALLOCATE(w(nmatp),rwork(7*nmatp))
  lwork=2*nmatp
  ALLOCATE(work(lwork))

  ! diagonalise the matrix
  CALL zhegvx(1,'V','I','U',nmatp,h,nmatp,o,nmatp,vl,vu,1,nstfv,evaltol,m,w, &
   evecfv,nmatmax,work,lwork,rwork,iwork,ifail,info)

  IF(info /= 0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(eveqnfvz): diagonalisation failed")')
    WRITE(*,'(" ZHEGVX RETURN ed INFO = ",I8)') info
    IF(info.gt.nmatp) THEN 
      i=info-nmatp
      WRITE(*,'(" The leading minor of the overlap matrix of order ",I8)') i
      WRITE(*,'("  is not positive definite")')
      WRITE(*,'(" Order of overlap matrix : ",I8)') nmatp
    ENDIF 
    WRITE(*,*)
    STOP 
  ENDIF 
  evalfv(1:nstfv)=w(1:nstfv)

  DEALLOCATE(iwork,ifail,w,rwork,work)

  CALL timesec(ts1)

  timefv=timefv+ts1-ts0
  RETURN 
END SUBROUTINE 

