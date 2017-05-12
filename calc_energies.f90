SUBROUTINE calc_energies( psi )
  USE m_states, ONLY : Nstates, Focc
  USE m_realspace, ONLY : dVol, Npoints
  USE m_PWGrid, ONLY : Ngwx
  USE m_hamiltonian, ONLY : V_ps_loc, Rhoe_R, V_Hartree
  USE m_energies
  IMPLICIT NONE 
  COMPLEX(8) :: psi(Ngwx,Nstates)
  !
  COMPLEX(8), ALLOCATABLE :: Kpsi(:)
  REAL(8), ALLOCATABLE :: epsxc(:)
  INTEGER :: ist
  COMPLEX(8) :: zdotc

  E_total   = 0.d0
  E_kinetic = 0.d0
  E_ps_loc  = 0.d0
  E_Hartree = 0.d0
  E_xc      = 0.d0

  ALLOCATE( Kpsi(Ngwx) )
  DO ist = 1,Nstates
    CALL op_K( 1, psi(:,ist), Kpsi(:) )
    E_kinetic = E_kinetic + Focc(ist) * real( zdotc( Ngwx, psi(:,ist),1, Kpsi(:),1 ), kind=8 )
  ENDDO 

  E_ps_loc = sum( Rhoe_R(:) * V_ps_loc(:) )*dVol

  E_Hartree = 0.5d0*sum( Rhoe_R(:) * V_Hartree(:) )*dVol

  ALLOCATE( epsxc(Npoints) )
  CALL excVWN( Npoints, Rhoe_R, epsxc )
  E_xc = sum( Rhoe_R(:) * epsxc(:) )*dVol

  E_total = E_kinetic + E_ps_loc + E_Hartree + E_xc + E_nn

  DEALLOCATE( epsxc )
  DEALLOCATE( Kpsi )
END SUBROUTINE 

