MODULE m_PsPot

  USE m_Ps_HGH, ONLY : Ps_HGH_Params_T
  IMPLICIT NONE 

  CHARACTER(64) :: PsPot_Dir = './HGH/'
  CHARACTER(64), ALLOCATABLE :: PsPot_FilePath(:)

  TYPE(Ps_HGH_Params_T), ALLOCATABLE :: Ps_HGH_Params(:)

  INTEGER :: NbetaNL
  REAL(8), ALLOCATABLE :: betaNL(:,:) ! (Npoints,NbetaNL)

  INTEGER :: NprojTotMax
  !REAL(8), ALLOCATABLE :: w_NL(:,:) ! Nspecies,0:3
  REAL(8), ALLOCATABLE :: w_NL(:)

END MODULE 


