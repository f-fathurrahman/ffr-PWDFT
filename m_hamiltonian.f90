MODULE m_hamiltonian

  IMPLICIT NONE

  REAL(8), ALLOCATABLE :: V_ps_loc(:)
  REAL(8), ALLOCATABLE :: V_Hartree(:)
  REAL(8), ALLOCATABLE :: V_xc(:)

  REAL(8), ALLOCATABLE :: Rhoe(:)

  REAL(8), ALLOCATABLE :: betaNL_psi(:,:,:) ! Natoms, NbetaNL, Nstates 

END MODULE 

