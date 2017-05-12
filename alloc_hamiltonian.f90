!! PURPOSE:
!!
!!   This subroutine allocates memory for global variables defined
!!   in module `m_hamiltonian`
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman

SUBROUTINE alloc_hamiltonian()
  USE m_hamiltonian, ONLY : V_ps_loc, V_Hartree, V_xc, Rhoe_R, betaNL_psi
  USE m_realspace, ONLY : Npoints
  USE m_PsPot, ONLY : NbetaNL
  USE m_atoms, ONLY : Natoms
  USE m_states, ONLY : Nstates
  IMPLICIT NONE

  ALLOCATE( V_ps_loc( Npoints ) )
  ALLOCATE( V_Hartree( Npoints ) )
  ALLOCATE( V_xc( Npoints ) )

  ALLOCATE( Rhoe_R( Npoints ) )

  ! XXX allocate here ???
  ALLOCATE( betaNL_psi(Natoms,Nstates,NbetaNL) )

  V_ps_loc(:) = 0.d0
  V_Hartree(:) = 0.d0
  V_xc(:) = 0.d0

  Rhoe_R(:) = 0.d0

END SUBROUTINE


