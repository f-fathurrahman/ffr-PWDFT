SUBROUTINE calc_energies( psi )
  USE m_states, ONLY : Nstates
  USE m_realspace, ONLY : dVol
  USE m_PWGrid, ONLY : Ngwx
  USE m_hamiltonian, ONLY : V_ps_loc
  USE m_energies
  IMPLICIT NONE 
  COMPLEX(8) :: psi(Ngwx,Nstates)
  !
  COMPLEX(8), ALLOCATABLE :: Kpsi(:)

  E_total   = 0.d0
  E_kinetic = 0.d0
  E_ps_loc  = 0.d0
  E_Hartree = 0.d0
  E_xc      = 0.d0

  ALLOCATE( Kpsi(Ngwx) )

  DO ist = 1,Nstates
    CALL op_K( 1, psi(:,ist), Kpsi(:) )
    E_kinetic = E_kinetic + Focc(ist) * zdotc( Ngwx, psi(:,ist),1, Kpsi(:),1 )
  ENDDO 

  E_ps_loc = sum( Rhoe(:) * V_ps_loc(:) )*dVol

END SUBROUTINE 

