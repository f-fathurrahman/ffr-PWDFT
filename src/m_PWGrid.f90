MODULE m_PWGrid

  IMPLICIT NONE 
  
  REAL(8) :: ecutwfc
  REAL(8) :: ecutrho

  INTEGER :: Ng
  REAL(8), ALLOCATABLE :: Gv(:,:)
  REAL(8), ALLOCATABLE :: Gv2(:)
  INTEGER, ALLOCATABLE :: idx_g2r(:)

  INTEGER :: Ngwx
  INTEGER, ALLOCATABLE :: idx_gw2g(:)
  INTEGER, ALLOCATABLE :: idx_gw2r(:)

END MODULE 

