SUBROUTINE KS_solve_SCF()
  USE m_constants, ONLY: PI
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol
  USE m_states, ONLY : Nstates, Focc, &
                       evecs => KS_evecs
  USE m_options, ONLY : ETHR_EVALS_LAST, ethr => ETHR_EVALS
  USE m_hamiltonian, ONLY : Rhoe
  USE m_energies, ONLY : Etot => E_total
  USE m_states, ONLY : Nelectrons
  USE m_options, ONLY : beta0, betamax, mixsdb, broydpm
  IMPLICIT NONE
  !
  INTEGER :: iterSCF
  REAL(8) :: Etot_old, dEtot
  REAL(8) :: dr2
  REAL(8), ALLOCATABLE :: Rhoe_old(:)

  REAL(8) :: integRho
  REAL(8), ALLOCATABLE :: beta_work(:), f_work(:)

  integRho = sum(Rhoe)*dVol
  WRITE(*,*)
  WRITE(*,'(1x,A,F18.10)') 'Initial guess: integRho = ', integRho

  ALLOCATE( Rhoe_old(Npoints) )

  beta0 = 0.1d0
  betamax = 1.d0
  ! Broyden parameters recommended by M. Meinert
  mixsdb = 5
  broydpm(1) = 0.4d0
  broydpm(2) = 0.15d0

  ALLOCATE( beta_work(Npoints) )
  ALLOCATE( f_work(Npoints) )

  Etot_old = 0.d0
  Rhoe_old(:) = Rhoe(:)

  !CALL mixadapt( 0, beta0, betamax, Npoints, Rhoe, Rhoe_old, beta_work, f_work, dr2 )
  
  dr2 = 1.d0
  DO iterSCF = 1, 100

    IF( iterSCF==1 ) THEN
      ethr = 1.d-1
    ELSE 
      IF( iterSCF == 2 ) ethr = 1.d-2
      !ethr = min( ethr, 1.d-2*dEtot / max(1.d0,Nelectrons) )
      ethr = ethr/5.d0
      ethr = max( ethr, ETHR_EVALS_LAST )
    ENDIF 

    CALL Sch_solve_diag()

    CALL calc_rhoe( Focc, evecs )

    CALL update_potentials()
    CALL calc_betaNL_psi( Nstates, evecs )
    CALL calc_energies( evecs ) ! update the potentials or not ?

    !CALL mixerifc( iterSCF, 1, Npoints, Rhoe, dr2, Npoints, Rhoe_old )
    CALL mixlinear( iterSCF, 0.1d0, Npoints, Rhoe, Rhoe_old, dr2 )

    !CALL mixadapt( iterSCF, beta0, betamax, Npoints, Rhoe, Rhoe_old, beta_work, f_work, dr2 )

    !IF( iterSCF > 2 ) THEN 
    !  dr2 = sqrt( ddot( Npoints, Rhoe(:)-Rhoe_old(:), 1, Rhoe(:)-Rhoe_old(:), 1 ) )
    !ENDIF

    integRho = sum(Rhoe)*dVol
    !WRITE(*,'(1x,A,F18.10)') 'After mix: integRho = ', integRho
    IF( abs(integRho - Nelectrons) > 1.0d-6 ) THEN
      WRITE(*,'(1x,A,ES18.10)') 'WARNING: diff after mix rho = ', abs(integRho-Nelectrons)
      WRITE(*,*) 'Rescaling Rho'
      Rhoe(:) = Nelectrons/integRho * Rhoe(:)
      integRho = sum(Rhoe)*dVol
      WRITE(*,'(1x,A,F18.10)') 'After rescaling: integRho = ', integRho
      WRITE(*,*)
    ENDIF 

!    CALL update_potentials()
!    CALL calc_betaNL_psi( Nstates, evecs )
!    CALL calc_energies( evecs ) ! update the potentials or not ?

    dEtot = abs(Etot - Etot_old)

    IF( dEtot < 1d-7) THEN 
      WRITE(*,*)
      WRITE(*,'(1x,A,I5,A)') 'SCF converged at ', iterSCF, ' iterations.'
      EXIT 
    ENDIF 

    !WRITE(*,*)
    WRITE(*,'(1x,A,I5,F18.10,2ES18.10)') 'SCF iter', iterSCF, Etot, dEtot, dr2

    Etot_old = Etot
    Rhoe_old(:) = Rhoe(:)
  ENDDO

  DEALLOCATE( Rhoe_old )
  DEALLOCATE( beta_work )
  DEALLOCATE( f_work )

END SUBROUTINE 

