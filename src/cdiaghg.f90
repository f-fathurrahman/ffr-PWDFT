!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! modified by Fadjar Fathurrahman
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE cdiaghg( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H hermitean matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... LAPACK version - uses both ZHEGV and ZHEGVX
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=8
  !
  INTEGER, INTENT(IN) :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculate
    ! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
    ! actually intent(in) but compilers don't know and complain
    ! matrix to be diagonalized
    ! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  COMPLEX(DP), INTENT(OUT) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  INTEGER                  :: lwork, nb, mm, info, i, j
    ! mm = number of calculated eigenvectors
  REAL(DP)                 :: abstol
  INTEGER,     ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    ALLOCATABLE :: rwork(:), sdiag(:), hdiag(:)
  COMPLEX(DP), ALLOCATABLE :: work(:)
    ! various work space
  LOGICAL                  :: all_eigenvalues
 ! REAL(DP), EXTERNAL       :: DLAMCH
  INTEGER,  EXTERNAL       :: ILAENV
    ! ILAENV returns optimal block size "nb"
  !
  !
  ! ... only the first processor diagonalizes the matrix
  !
     !
     ! ... save the diagonal of input S (it will be overwritten)
     !
     ALLOCATE( sdiag( n ) )
     DO i = 1, n
        sdiag(i) = DBLE( s(i,i) )
     END DO
     !
     all_eigenvalues = ( m == n )
     !
     ! ... check for optimal block size
     !
     nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
     !
     IF ( nb < 1 .OR. nb >= n) THEN
        !
        lwork = 2*n
        !
     ELSE
        !
        lwork = ( nb + 1 )*n
        !
     END IF
     !
     ALLOCATE( work( lwork ) )
     !
     IF ( all_eigenvalues ) THEN
        !
        ALLOCATE( rwork( 3*n - 2 ) )
        !
        ! ... calculate all eigenvalues (overwritten to v)
        !
        v(:,:) = h(:,:)
        !
        CALL ZHEGV( 1, 'V', 'U', n, v, ldh, &
                    s, ldh, e, work, lwork, rwork, info )
        !
     ELSE
        !
        ALLOCATE( rwork( 7*n ) )
        !
        ! ... save the diagonal of input H (it will be overwritten)
        !
        ALLOCATE( hdiag( n ) )
        DO i = 1, n
           hdiag(i) = DBLE( h(i,i) )
        END DO
        !
        ALLOCATE( iwork( 5*n ) )
        ALLOCATE( ifail( n ) )
        !
        ! ... calculate only m lowest eigenvalues
        !
        abstol = 0.D0
       ! abstol = 2.D0*DLAMCH( 'S' )
        !
        ! ... the following commented lines calculate optimal lwork
        !
        !lwork = -1
        !
        !CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
        !             0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
        !             work, lwork, rwork, iwork, ifail, info )
        !
        !lwork = INT( work(1) ) + 1
        !
        !IF( lwork > SIZE( work ) ) THEN
        !   DEALLOCATE( work )
        !   ALLOCATE( work( lwork ) )
        !END IF
        !
        CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                     work, lwork, rwork, iwork, ifail, info )
        !
        DEALLOCATE( ifail )
        DEALLOCATE( iwork )
        !
        ! ... restore input H matrix from saved diagonal and lower triangle
        !
        DO i = 1, n
           h(i,i) = CMPLX( hdiag(i), 0.0_DP ,kind=DP)
           DO j = i + 1, n
              h(i,j) = CONJG( h(j,i) )
           END DO
           DO j = n + 1, ldh
              h(j,i) = ( 0.0_DP, 0.0_DP )
           END DO
        END DO
        !
        DEALLOCATE( hdiag )
        !
     END IF
     !
     !
     DEALLOCATE( rwork )
     DEALLOCATE( work )
     !
     IF ( info > n ) THEN
       WRITE(*,*) 'ERROR in cdiaghg: S matrix is not positive definite'
       STOP 
     ELSE IF ( info > 0 ) THEN
       WRITE(*,*) 'ERROR in cdiaghg: eigenvectors failed to converge'
       STOP 
     ELSE IF ( info < 0 ) THEN
       WRITE(*,*) 'ERROR in cdiaghg: incorrect call to ZHEGV*'
       STOP 
     END IF
     !
     ! ... restore input S matrix from saved diagonal and lower triangle
     !
     DO i = 1, n
        s(i,i) = CMPLX( sdiag(i), 0.0_DP ,kind=DP)
        DO j = i + 1, n
           s(i,j) = CONJG( s(j,i) )
        END DO
        DO j = n + 1, ldh
           s(j,i) = ( 0.0_DP, 0.0_DP )
        END DO
     END DO
     !
     DEALLOCATE( sdiag )
     !
  !
  RETURN
  !
END SUBROUTINE cdiaghg
!
