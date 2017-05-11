MODULE m_options

  IMPLICIT NONE 

  ! Options for controlling how beta is calculated
  INTEGER :: CG_BETA = 2
  ! 1 => Fletcher-Reeves
  ! 2 => Polak-Ribiere
  ! 3 => Hestenes-Stiefel
  ! 4 => Dai-Yuan

  ! whether free nabla2 after constructing matrix or not
  LOGICAL :: FREE_NABLA2 = .FALSE.

  REAL(8) :: DIAG_DAVIDSON_QE_ETHR = 1.0d-5

  ! type of mixing to use for the potential
  integer :: mixtype
  ! mixing type description
  character(256) :: mixdescr
  ! adaptive mixing parameter
  REAL(8) :: beta0
  REAL(8) :: betamax
  ! subspace dimension for Broyden mixing
  INTEGER :: mixsdb
  ! Broyden mixing parameters alpha and w0
  REAL(8) :: broydpm(2)


END MODULE 

