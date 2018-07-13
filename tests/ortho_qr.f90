!-------------------------------------------
SUBROUTINE z_ortho_qr(X, nbasis, nstates)
!-------------------------------------------
  IMPLICIT NONE
  ! Argument
  INTEGER :: nbasis, nstates
  COMPLEX(8) :: X(nbasis,nstates)
  ! Local
  INTEGER :: m, n, k, lwork, info
  COMPLEX(8), ALLOCATABLE :: tau(:), work(:)

  m = nbasis
  n = nstates
  k = nstates
  lwork = n

  ALLOCATE(tau(k))
  ALLOCATE(work(lwork))

  ! Compute QR factorization
  CALL zgeqrf(m,n,X,m,tau,work,lwork,info)
  IF(info /= 0) THEN 
    WRITE (*,'(1x,A,I4)') 'ERROR calling ZGEQRF : info = ', info
    RETURN 
  ENDIF

  ! Construct Q explicitly
  CALL zungqr(m,k,k,X,m,tau,work,lwork,info)
  if(info /= 0) then
    WRITE(*,'(1x,A,I4)') 'ERROR calling ZUNGQR : info = ', info
    RETURN 
  ENDIF

  DEALLOCATE(tau)
  DEALLOCATE(work)
END  SUBROUTINE 

