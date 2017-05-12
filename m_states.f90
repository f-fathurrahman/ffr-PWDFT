MODULE m_states

  IMPLICIT NONE 
  
  INTEGER :: Nstates
  REAL(8) :: Nelectrons

  REAL(8), ALLOCATABLE :: KS_evals(:)
  REAL(8), ALLOCATABLE :: KS_evecs(:,:)

  REAL(8), ALLOCATABLE :: Focc(:)

END MODULE 

