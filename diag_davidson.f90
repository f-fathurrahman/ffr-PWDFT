SUBROUTINE diag_davidson( evals, V, TOLERANCE )
  USE m_states, only : Nstates
  USE m_PWGrid, only : Ngwx
  IMPLICIT NONE
  COMPLEX(8), INTENT(inout) :: V(Ngwx,Nstates)
  REAL(8), INTENT(inout) :: EVALS(Nstates)
  ! Local variable
  REAL(8) :: RNORM
  INTEGER :: ist, istep, MAX_DIR, NCONV
  LOGICAL :: IS_CONVERGED
  REAL(8) :: MACHINE_ZERO, TOLERANCE
  REAL(8), ALLOCATABLE :: RES_TOL(:), RES_NORM(:), EVALS_RED(:)
  COMPLEX(8) :: Z_ZERO, Z_ONE
  COMPLEX(8), ALLOCATABLE :: CMAT(:,:), H_MAT(:,:), O_MAT(:,:), EVEC(:,:)
  COMPLEX(8), ALLOCATABLE :: HV(:,:), R(:,:), HR(:,:), XTEMP(:,:)
  ! BLAS function
  COMPLEX(8) :: ZDOTC

  Z_ZERO = CMPLX(0.D0,0.D0)
  Z_ONE = CMPLX(1.D0,0.D0)

  ALLOCATE(RES_TOL(Nstates))
  ALLOCATE(RES_NORM(Nstates))
  ALLOCATE(CMAT(Nstates,Nstates))
  ALLOCATE(H_MAT(2*Nstates,2*Nstates))
  ALLOCATE(O_MAT(2*Nstates,2*Nstates))
  ALLOCATE(EVEC(2*Nstates,2*Nstates))
  ALLOCATE(EVALS_RED(2*Nstates))
  ALLOCATE(HV(Ngwx,Nstates))
  ALLOCATE(R(Ngwx,Nstates))
  ALLOCATE(HR(Ngwx,Nstates))
  ALLOCATE(XTEMP(Ngwx,Nstates))

  WRITE(*,*) 'diag_davidson pass 34'

!  V(:,:) = Z_ZERO
!  DO ist=1,Nstates
!    V(ist,ist) = Z_ONE;
!  ENDDO

  WRITE(*,*) 'V, HV = ', sum(V), sum(HV)
  ! Apply Hamiltonian
  call op_H( Nstates, V, HV )

  WRITE(*,*) 'diag_davidson pass 44'

  ! Calculate Rayleigh quotient
  WRITE(*,*) 'shape(V) = ', shape(V)
  WRITE(*,*) 'shape(HV) = ', shape(HV)
  WRITE(*,*) 'shape(evals) = ', shape(evals)
  WRITE(*,*) 'V, HV = ', sum(V), sum(HV)
  !WRITE(*,*) 'V = ', V
  DO ist=1,Nstates
    WRITE(*,*) 'ist = ', ist
    EVALS(ist) = REAL( ZDOTC(Ngwx, V(:,ist),1, HV(:,ist),1),kind=8 )
    !WRITE(*,*) 'evals = ', evals(ist)
  ENDDO
!  stop

  WRITE(*,*) 'diag_davidson pass 52'

  ! Calculate matrix of residual vector
  DO ist=1,Nstates
    R(:,ist) = EVALS(ist)*V(:,ist) - HV(:,ist)
    RES_TOL(ist) = SQRT( real(ZDOTC(Ngwx, R(:,ist),1, R(:,ist),1),kind=8) )
  ENDDO

  WRITE(*,*) 'diag_davidson pass 60'

  istep = 1
  IS_CONVERGED = .FALSE.
  MAX_DIR = 100
  MACHINE_ZERO = 2.220446049250313D-16
  RNORM = 1.D0

  DO WHILE ( (istep <= MAX_DIR) .AND. (.NOT.IS_CONVERGED) )
    WRITE(*,'(2I8,F18.10)') istep, NCONV, RNORM
    RES_NORM = 1.D0

    ! WHERE(MACHINE_ZERO < RES_TOL) RES_NORM = 1.0_DP/RES_TOL
    !WRITE(*,*) 'RES_NORM:'
    DO ist = 1,Nstates
      IF(MACHINE_ZERO < RES_TOL(ist)) RES_NORM(ist) = 1.D0/RES_TOL(ist)
      !WRITE(*,*) RES_NORM(ist)
    END DO

    ! Scale the residual vectors
    DO ist=1,Nstates
      R(:,ist) = RES_NORM(ist)*R(:,ist)
    ENDDO

    ! Apply preconditioner
    call prec_Gv2_inplace(Nstates,R)


! Construct the reduced hamiltonian. The reduced hamiltonian has dimensions
!  2nb x 2nb and is constructed by filling in four nb x nb blocks one at a time:
! __
!|  |
!| <v|H|v>   <v|H|r>  |
!    h_mat = |  |
!| *******   <r|H|r>  |
!|__|

    call op_H(Nstates,R,HR)

    IF(istep == 1) THEN
      CALL ZGEMM('C','N',Nstates,Nstates,Ngwx, Z_ONE,V,Ngwx, HV,Ngwx, Z_ZERO,CMAT,Nstates)
      H_MAT(1:Nstates,1:Nstates) = CMAT
    ELSE
      H_MAT(1:Nstates,1:Nstates) = Z_ZERO
      DO ist = 1,Nstates
        H_MAT(ist,ist) = CMPLX(EVALS(ist),0.D0)
      ENDDO
    ENDIF
    ! <v|H|r> --> cmat
    CALL ZGEMM('C','N',Nstates,Nstates,Ngwx, Z_ONE,V,Ngwx, HR,Ngwx, Z_ZERO,CMAT,Nstates)
    H_MAT(1:Nstates,Nstates+1:2*Nstates) = CMAT
    CALL ZGEMM('C','N',Nstates,Nstates,Ngwx, Z_ONE,HR,Ngwx, V,Ngwx, Z_ZERO,CMAT,Nstates)
    H_MAT(Nstates+1:2*Nstates,1:Nstates) = CMAT
    ! <r|H|r> --> cmat
    CALL ZGEMM('C','N',Nstates,Nstates,Ngwx, Z_ONE,R,Ngwx, HR,Ngwx, Z_ZERO,CMAT,Nstates)
    H_MAT(Nstates+1:2*Nstates,Nstates+1:2*Nstates) = CMAT

! Construct the reduced overlap matrix which has dimenstions 2nb x 2nb
!   and is constructed by filling in four nb x nb blocks one at a time:
! _   _
!|     |
!|  <v|v>   <v|r>  |
!    o_mat = |     |
!|  *****   <r|r>  |
!|_   _|

    O_MAT(1:Nstates,1:Nstates) = Z_ZERO
    DO ist = 1,Nstates
      O_MAT(ist,ist) = Z_ONE
    END DO
    ! <v|r> --> cmat
    CALL ZGEMM('C','N',Nstates,Nstates,Ngwx, Z_ONE,V,Ngwx, R,Ngwx, Z_ZERO,CMAT,Nstates)
    O_MAT(1:Nstates,Nstates+1:2*Nstates) = CMAT
    CALL ZGEMM('C','N',Nstates,Nstates,Ngwx, Z_ONE,R,Ngwx, V,Ngwx, Z_ZERO,CMAT,Nstates)
    O_MAT(Nstates+1:2*Nstates,1:Nstates) = CMAT
    ! <r|r> --> cmat
    CALL ZGEMM('C','N',Nstates,Nstates,Ngwx, Z_ONE,R,Ngwx, R,Ngwx, Z_ZERO,CMAT,Nstates)
    O_MAT(Nstates+1:2*Nstates,Nstates+1:2*Nstates) = CMAT

    !CALL EIG_ZHEGV(H_MAT,O_MAT,EVALS_RED)
    CALL EIG_ZHEGV_EVAL_F90(H_MAT,2*Nstates, O_MAT,2*Nstates, EVALS_RED, EVEC,2*Nstates, 2*Nstates)
    !CALL EIG_ZHEGV(H_MAT,O_MAT,EVALS_RED,2*Nstates)

    !WRITE(*,*) 'REDUCED EVALS = '
    !DO I=1,2*Nstates
    !  WRITE(*,*) EVALS_RED(I)
    !ENDDO

    EVALS = EVALS_RED(1:Nstates)
    CMAT = EVEC(1:Nstates,1:Nstates)
    ! V*CMAT --> V
    CALL ZGEMM('N','N',Ngwx,Nstates,Nstates, Z_ONE,V,Ngwx, CMAT,Nstates, Z_ZERO,XTEMP,Ngwx)
    V = XTEMP
    ! HV = HV*CMAT
    CALL ZGEMM('N','N',Ngwx,Nstates,Nstates, Z_ONE,HV,Ngwx, CMAT,Nstates, Z_ZERO,XTEMP,Ngwx)
    HV = XTEMP
    !
    CMAT = EVEC(Nstates+1:2*Nstates,1:Nstates)
    ! V = V + R*CMAT
    CALL ZGEMM('N','N',Ngwx,Nstates,Nstates, Z_ONE,R,Ngwx, CMAT,Nstates, Z_ONE,V,Ngwx)
    ! HV = HV + HR*CMAT
    CALL ZGEMM('N','N',Ngwx,Nstates,Nstates, Z_ONE,HR,Ngwx, CMAT,Nstates, Z_ONE,HV,Ngwx)

    ! Calculate matrix of residual vector
    DO ist=1,Nstates
      R(:,ist) = EVALS(ist)*V(:,ist) - HV(:,ist)
      RES_TOL(ist) = SQRT( real(ZDOTC(Ngwx, R(:,ist),1, R(:,ist),1),kind=8) )
      !WRITE(*,'(1X,I5,2F18.10)') ist, EVALS(ist), RES_TOL(ist)
    ENDDO

    IS_CONVERGED = .TRUE.
    DO ist = 1,Nstates
      IS_CONVERGED = (IS_CONVERGED .AND. (RES_TOL(ist) < TOLERANCE) )
    END DO
    istep = istep + 1
    RNORM = SUM(RES_TOL)/REAL(Nstates,8)
  END DO

  !rnorm = sum(res_tol)/real(nstates,8)

  DO ist=1,Nstates
    WRITE(*,'(1X,I5,2F18.10)') ist, EVALS(ist), RES_TOL(ist)
  ENDDO
  WRITE(*,*) 'END OF DAVIDSON ITERATION: RNORM = ', RNORM

  DEALLOCATE(RES_TOL)
  DEALLOCATE(RES_NORM)
  DEALLOCATE(CMAT)
  DEALLOCATE(H_MAT)
  DEALLOCATE(O_MAT)
  DEALLOCATE(EVEC)
  DEALLOCATE(EVALS_RED)
  DEALLOCATE(HV)
  DEALLOCATE(R)
  DEALLOCATE(HR)
  DEALLOCATE(XTEMP)
END SUBROUTINE


!---------------------------------------------------
SUBROUTINE EIG_ZHEGV_EVAL_F90(A,LDA, B,LDB, EVAL, EVEC,LDE, N)
!---------------------------------------------------
  IMPLICIT NONE
  ! ARGUMENTS
  INTEGER :: LDA,LDB,LDE,N
  COMPLEX(8) :: A(LDA,N)
  COMPLEX(8) :: B(LDB,N)
  COMPLEX(8) :: EVEC(LDE,N)
  REAL(8) :: EVAL(N)
  ! LOCAL
  REAL(8), PARAMETER :: SMALL=2.220446049250313D-16
  COMPLEX(8), PARAMETER :: Z_ZERO=(0.D0,0.D0), Z_ONE=(1.D0,0.D0)
  INTEGER :: LWORK, LRWORK, LIWORK, INFO, I, NN
  COMPLEX(8), ALLOCATABLE :: WORK(:)
  REAL(8), ALLOCATABLE :: RWORK(:)
  INTEGER, ALLOCATABLE :: IWORK(:)
  REAL(8) :: SCAL

  LWORK = N*N + 2*N
  LRWORK = 2*N*N + 5*N + 1
  LIWORK = 5*N + 3

  ALLOCATE(WORK(LWORK))
  ALLOCATE(RWORK(LRWORK))
  ALLOCATE(IWORK(LIWORK))

  ! DIAGONALIZE B
  CALL ZHEEVD('V','U',N, B,LDB, EVAL, WORK,LWORK, RWORK,LRWORK, IWORK,LIWORK, INFO)
  IF(INFO /= 0) THEN
    WRITE(*,'(1X,A,I4)') 'ERROR CALLING ZHEEVD IN EIG_ZHEGV_F90 : INFO = ', INFO
    STOP
  ENDIF

  NN = 0
  DO I=1,N
    IF( EVAL(I) > SMALL ) THEN
      NN = NN + 1
      SCAL = 1.D0/SQRT(EVAL(I))
      CALL ZDSCAL(N, SCAL, B(:,I),1)
    ENDIF
  ENDDO
  IF(NN < N) THEN
    WRITE(*,'(1X,A,I4)') 'WARNING: NUMBER OF LINEARLY INDEPENDENT VECTORS = ', NN
    WRITE(*,'(1X,A,I4)') '       WHILE SIZE OF THE PROBLEM = ', N
  ENDIF

  ! TRANSFORM A:
  ! A <-- EVEC(B)* A EVEC(B)
  CALL ZGEMM('N','N', N,N,N, Z_ONE,A,LDA, B,LDB, Z_ZERO,EVEC,LDE)
  CALL ZGEMM('C','N', N,N,N, Z_ONE,B,LDB, EVEC,LDE, Z_ZERO,A,LDA)

  ! DIAGONALIZE TRANSFORMED A
  CALL ZHEEVD('V','U',N, A,LDA, EVAL, WORK,LWORK, RWORK,LRWORK, IWORK,LIWORK, INFO)
  IF(INFO /= 0) THEN
    WRITE(*,'(1X,A,I4)') 'ERROR CALLING ZHEEVD IN EIG_ZHEEVD_F90 : INFO = ', INFO
    STOP
  ENDIF

  ! BACK TRANSFORM EIGENVECTORS
  CALL ZGEMM('N','N', N,N,N, Z_ONE,B,LDB, A,LDA, Z_ZERO,EVEC,LDE)

  DEALLOCATE(WORK)
  DEALLOCATE(RWORK)
  DEALLOCATE(IWORK)
END SUBROUTINE
