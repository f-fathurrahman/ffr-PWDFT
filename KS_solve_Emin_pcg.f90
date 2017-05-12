!! PURPOSE:
!!
!!   This subroutine solves Kohn-Sham equations by minimizing total energy
!!   functional using conjugate gradient algorithm.
!!   The algorithm is based on T.A. Arias notes.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! MODIFIES:
!! 
!!   Global variables `KS_evecs` and `E_total`
!!
!! NOTES:
!!
!!   ILU0 preconditioner from SPARSKIT is used as preconditioner.
!!

SUBROUTINE KS_solve_Emin_pcg( alpha_t, NiterMax, restart )

  USE m_constants, ONLY : ZZERO
  USE m_PWGrid, ONLY : Ngwx
  USE m_states, ONLY : Nstates, &
                       Focc, &
                       v => KS_evecs
  USE m_energies, ONLY : Etot => E_total 
  USE m_options, ONLY : CG_BETA

  IMPLICIT NONE
  !
  INTEGER :: NiterMax
  REAL(8) :: alpha_t  ! step size
  LOGICAL :: restart
  COMPLEX(8), ALLOCATABLE :: g(:,:), g_old(:,:), g_t(:,:)
  COMPLEX(8), ALLOCATABLE :: d(:,:), d_old(:,:)
  COMPLEX(8), ALLOCATABLE :: Kg(:,:), Kg_old(:,:) ! preconditioned
  COMPLEX(8), ALLOCATABLE :: tv(:,:)
  REAL(8) :: alpha, beta, denum, Etot_old
  !
  INTEGER :: iter

  CALL info_KS_solve_Emin_pcg( alpha_t, NiterMax, restart )

  ALLOCATE( g(Ngwx,Nstates) )
  ALLOCATE( g_old(Ngwx,Nstates) )
  ALLOCATE( g_t(Ngwx,Nstates) )
  ALLOCATE( d(Ngwx,Nstates) )
  ALLOCATE( d_old(Ngwx,Nstates) )

  ALLOCATE( Kg(Ngwx,Nstates) )
  ALLOCATE( Kg_old(Ngwx,Nstates) )

  ALLOCATE( tv(Ngwx,Nstates) )


  ! Read starting eigenvectors from file
  IF( restart ) THEN
    READ(112) v   ! FIXME Need to use file name
  ENDIF

  CALL calc_Rhoe_R( Focc, v )
!  CALL update_potentials()
!  CALL calc_betaNL_psi( Nstates, v )
  CALL calc_energies( v )
  !CALL info_energies()

  Etot_old = Etot

  alpha = 0.d0
  beta  = 0.d0

  g(:,:)      = ZZERO
  g_t(:,:)    = ZZERO
  d(:,:)      = ZZERO
  d_old(:,:)  = ZZERO
  Kg(:,:)     = ZZERO
  Kg_old(:,:) = ZZERO

  DO iter = 1, NiterMax
    !
    ! Evaluate gradient at current trial vectors
    CALL calc_grad( Nstates, v, g )
    ! Precondition
    CALL prec_Gv2( Nstates, g(:,:), Kg(:,:) )
    !
    ! set search direction
    IF( iter /= 1 ) THEN
      SELECT CASE ( CG_BETA )
      CASE(1)
        ! Fletcher-Reeves
        beta = real( sum( conjg(g) * Kg ) / sum( conjg(g_old) * Kg_old ), kind=8 )
      CASE(2)
        ! Polak-Ribiere
        beta = real( sum( conjg(g-g_old)*Kg ) / sum( conjg(g_old) * Kg_old ), kind=8 )
      CASE(3)
        ! Hestenes-Stiefel
        beta = real( sum( conjg(g-g_old)*Kg ) / sum( conjg(g-g_old)*d_old ), kind=8 )
      CASE(4)
        ! Dai-Yuan
        beta = real( sum( conjg(g) * Kg ) / sum( conjg(g-g_old)*d_old ), kind=8 )
      END SELECT 
    ENDIF
    IF( beta < 0 ) THEN 
      WRITE(*,'(1x,A,F18.10,A)') 'beta is smaller than zero: ', beta, ': setting it to zero'
    ENDIF 
    beta = max( 0.d0, beta )
    d(:,:) = -Kg(:,:) + beta*d_old(:,:)
    !
    ! Evaluate gradient at trial step
    tv(:,:) = v(:,:) + alpha_t * d(:,:)
    CALL z_ortho_gram_schmidt( tv, Ngwx, Ngwx, Nstates )

    CALL calc_Rhoe_R( Focc, tv )
!    CALL update_potentials()  ! Now global vars on m_hamiltonian are changed
!    CALL calc_betaNL_psi( Nstates, tv )
    CALL calc_grad( Nstates, tv, g_t )
    !
    ! Compute estimate of best step and update current trial vectors
    denum = real( sum( conjg(g - g_t) * d ), kind=8 )
    IF( denum /= 0.d0 ) THEN  ! FIXME: use abs ?
      alpha = abs( alpha_t * real(sum(conjg(g)*d),kind=8)/denum )
    ELSE 
      alpha = 0.d0
    ENDIF
    !WRITE(*,*) 'iter, alpha, beta', iter, alpha, beta

    v(:,:) = v(:,:) + alpha * d(:,:)
    CALL z_ortho_gram_schmidt( v, Ngwx, Ngwx, Nstates )

    CALL calc_Rhoe_R( Focc, v )
!    CALL update_potentials()
!    CALL calc_betaNL_psi( Nstates, v )
    CALL calc_energies( v )
    !
    WRITE(*,'(1x,I5,F18.10,ES18.10)') iter, Etot, Etot_old-Etot
    !
    IF( abs(Etot - Etot_old) < 1.d-7 ) THEN
      WRITE(*,*) 'KS_solve_Emin_pcg converged in iter', iter
      EXIT
    ENDIF
    !
    Etot_old = Etot
    g_old(:,:) = g(:,:)
    d_old(:,:) = d(:,:)
    Kg_old(:,:) = Kg(:,:)
  ENDDO

  WRITE(111) v

  DEALLOCATE( g, g_old, g_t, d, d_old, tv, Kg, Kg_old )
END SUBROUTINE


SUBROUTINE info_KS_solve_Emin_pcg( alpha_t, NiterMax, restart )
  USE m_options, ONLY : CG_BETA
  USE m_PWGrid, ONLY : Ngwx
  USE m_states, ONLY : Nstates
  IMPLICIT NONE 
  INTEGER :: NiterMax
  REAL(8) :: alpha_t
  LOGICAL :: restart
  !
  REAL(8) :: memGB

  memGB = Ngwx*Nstates*8d0 * 8d0 / (1024d0*1024d0*1024.d0)

  WRITE(*,*)
  WRITE(*,*) 'Minimizing KS total energy functional using PCG algorithm'
  WRITE(*,*) '---------------------------------------------------------'
  WRITE(*,*)
  WRITE(*,*) 'NiterMax = ', NiterMax
  WRITE(*,*) 'alpha_t  = ', alpha_t
  WRITE(*,*) 'restart  = ', restart
  WRITE(*,*)
  IF( CG_BETA == 1 ) THEN
    WRITE(*,*) 'Using Fletcher-Reeves formula'
  ELSEIF( CG_BETA == 2 ) THEN 
    WRITE(*,*) 'Using Polak-Ribiere formula'
  ELSEIF( CG_BETA == 3 ) THEN 
    WRITE(*,*) 'Using Hestenes-Stiefel formula'
  ELSEIF( CG_BETA == 4 ) THEN 
    WRITE(*,*) 'Using Dai-Yuan formula'
  ELSE 
    WRITE(*,*) 'XXXXX WARNING: Unknow CG_BETA: ', CG_BETA
  ENDIF 
  WRITE(*,*)
  WRITE(*,*) 'KS_solve_Emin_pcg: memGB = ', memGB
  WRITE(*,*)

END SUBROUTINE 
