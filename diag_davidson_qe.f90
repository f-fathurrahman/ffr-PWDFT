!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! modified by Fadjar Fathurrahman from cegterg.f90 of QE-6.1 distribution.
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!
!----------------------------------------------------------------------------
SUBROUTINE diag_davidson_qe( npw, npwx, nvec, nvecx, evc, ethr, &
                             e, btype, notcnv, dav_iter )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector
  !
  IMPLICIT NONE
  integer, parameter :: DP=8
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set :
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! umber of spin polarizations
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx,nvec)
    !  evc contains the  refined estimates of the eigenvectors
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence :
    !   root improvement is stopped, when two consecutive estimates of the root
    !   differ by less than ethr.
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer number of iterations performed
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 100
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! adapted npw and npwx
    ! do-loop counters
  INTEGER :: ierr
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  REAL(DP), ALLOCATABLE :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr
    ! threshold for empty bands
  !
  REAL(DP), EXTERNAL :: ddot
  !
  ! EXTERNAL  h_psi,    s_psi,    g_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  if( nvec > nvecx/2 ) then
    write(*,*) 'Error in diag_davidson_qe: nvecx is too small:', nvecx
    stop
  endif
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  !
  kdim = npw
  kdmx = npwx
  !
  ALLOCATE( psi( npwx, nvecx ) )
  ALLOCATE( hpsi( npwx, nvecx ) )
  !
  ALLOCATE( sc( nvecx, nvecx ) )
  ALLOCATE( hc( nvecx, nvecx ) )
  ALLOCATE( vc( nvecx, nvecx ) )
  ALLOCATE( ew( nvecx ) )
  ALLOCATE( conv( nvec ) )
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  hpsi = ZERO
  psi  = ZERO
  psi(:,1:nvec) = evc(:,1:nvec)
  !
  ! ... hpsi contains h times the basis vectors
  !
  call op_H( nvec, psi, hpsi )
  !
  ! ... hc contains the projection of the hamiltonian onto the reduced
  ! ... space vc contains the eigenvectors of hc
  !
  hc(:,:) = ZERO
  sc(:,:) = ZERO
  vc(:,:) = ZERO
  !
  CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
              psi, kdmx, hpsi, kdmx, ZERO, hc, nvecx )
  CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
                 psi, kdmx, psi, kdmx, ZERO, sc, nvecx )
  !
  ! ... diagonalize the reduced hamiltonian
  !
  CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
  !
  e(1:nvec) = ew(1:nvec)
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter
     !
     np = 0
     !
     DO n = 1, nvec
        !
        IF ( .NOT. conv(n) ) THEN
           !
           ! ... this root not yet converged ...
           !
           np = np + 1
           !
           ! ... reorder eigenvectors so that coefficients for unconverged
           ! ... roots come first. This allows to use quick matrix-matrix
           ! ... multiplications to set a new basis vector (see below)
           !
           IF ( np /= n ) vc(:,np) = vc(:,n)
           !
           ! ... for use in g_psi
           !
           ew(nbase+np) = e(n)
           !
        END IF
        !
     END DO
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, psi, &
                kdmx, vc, nvecx, ZERO, psi(1,nb1), kdmx )
     !
     DO np = 1, notcnv
        psi(:,nbase+np) = - ew(nbase+np)*psi(:,nbase+np)
     END DO
     !
     CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, hpsi, &
                 kdmx, vc, nvecx, ONE, psi(1,nb1), kdmx )
     !
     ! ... approximate inverse iteration
     !
     !CALL g_psi( npwx, npw, notcnv, npol, psi(1,1,nb1), ew(nb1) )
     CALL prec_Gv2_inplace( notcnv, psi(:,nb1) )
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
     ! ... order to improve numerical stability of subspace diagonalization
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     DO n = 1, notcnv
       !
       nbn = nbase + n
       !
       ew(n) = ddot( 2*npw, psi(1,nbn), 1, psi(1,nbn), 1 )
       !
     END DO
     !
     DO n = 1, notcnv
        psi(:,nbase+n) = psi(:,nbase+n) / SQRT( ew(n) )
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     !CALL h_psi( npwx, npw, notcnv, psi(1,1,nb1), hpsi(1,1,nb1) )
     call op_H( notcnv, psi(:,nb1), hpsi(:,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
                 kdmx, hpsi(1,nb1), kdmx, ZERO, hc(1,nb1), nvecx )
     !
     !
     CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
                 kdmx, psi(1,nb1), kdmx, ZERO, sc(1,nb1), nvecx )
     !
     nbase = nbase + notcnv
     !
     DO n = 1, nbase
        !
        ! ... the diagonal of hc and sc must be strictly real
        !
        hc(n,n) = CMPLX( REAL( hc(n,n) ), 0.D0 ,kind=DP)
        sc(n,n) = CMPLX( REAL( sc(n,n) ), 0.D0 ,kind=DP)
        !
        DO m = n + 1, nbase
           !
           hc(m,n) = CONJG( hc(n,m) )
           sc(m,n) = CONJG( sc(n,m) )
           !
        END DO
        !
     END DO
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
     !
     ! ... test for convergence
     !
     WHERE( btype(1:nvec) == 1 )
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
        !
     END WHERE
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     e(1:nvec) = ew(1:nvec)
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, &
                    psi, kdmx, vc, nvecx, ZERO, evc, kdmx )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
           !!!     &   " eigenvalues not converged")' ) notcnv
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        psi(:,1:nvec) = evc(:,1:nvec)
        !
        CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, hpsi, &
                    kdmx, vc, nvecx, ZERO, psi(1,nvec+1), kdmx )
        !
        hpsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        hc(:,1:nbase) = ZERO
        sc(:,1:nbase) = ZERO
        vc(:,1:nbase) = ZERO
        !
        DO n = 1, nbase
           !
           hc(n,n) = CMPLX( e(n), 0.0_DP ,kind=DP)
           !
           sc(n,n) = ONE
           vc(n,n) = ONE
           !
        END DO
        !
     END IF
     !
  END DO iterate
  !
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vc )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )
  !
  RETURN
  !
END SUBROUTINE
