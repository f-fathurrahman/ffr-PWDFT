!! PURPOSE:
!!
!!   This subroutine solves Schrodinger equation using iterative
!!   diagonalization.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! MODIFIES:
!!
!!   Global variables `KS_evecs` and `KS_evals`.
!!
!! IMPORTANT
!!
!!   `KS_evecs` should be initialized outside this subroutine.
!!   In the subsequent calls, `KS_evecs` is used as initial guesses.

SUBROUTINE Sch_solve_diag()

  use m_realspace, only : dVol
  USE m_states, ONLY : Nstates, Focc, &
                       evecs => KS_evecs, &
                       evals => KS_evals
  USE m_options, ONLY : ethr => DIAG_DAVIDSON_QE_ETHR
  USE m_PWGrid, ONLY : Ngwx
  IMPLICIT NONE
  INTEGER, ALLOCATABLE :: btype(:)
  INTEGER :: dav_iter
  INTEGER :: notcnv
  INTEGER :: ist

  ALLOCATE( btype(Nstates) )
  btype(:) = 1 ! assume all bands are occupied
  DO ist = 1,Nstates
    IF( Focc(ist) <= 1d-13 ) btype(ist) = 0
  ENDDO

  !WRITE(*,*)
  !WRITE(*,*) 'Solving Schrodinger equation with Davidson iterative diagonalization'
  !WRITE(*,*)

  CALL diag_davidson_qe( Ngwx, Ngwx, Nstates, 3*Nstates, evecs, ethr, &
                         evals, btype, notcnv, dav_iter )
  
  
  !CALL diag_davidson( evals, evecs, ethr )

  !CALL diag_lobpcg( Nstates, evals, evecs )

  !WRITE(*,'(1x,A,ES18.10,A,I4)') 'Davidson_QE: ethr = ', ethr, ' dav_iter = ', dav_iter

  WRITE(*,*)
  WRITE(*,*) 'Eigenvalues:'
  WRITE(*,*)
  DO ist = 1, Nstates
    WRITE(*,'(1x,I4,F18.10)') ist, evals(ist)
  ENDDO

  ! normalize evecs properly
  evecs(:,:) = evecs(:,:)/sqrt(dVol)

  DEALLOCATE( btype )
END SUBROUTINE
