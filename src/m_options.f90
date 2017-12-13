MODULE m_options

  IMPLICIT NONE 

  ! Options for controlling how beta is calculated
  INTEGER :: CG_BETA = 2
  ! 1 => Fletcher-Reeves
  ! 2 => Polak-Ribiere
  ! 3 => Hestenes-Stiefel
  ! 4 => Dai-Yuan

  REAL(8) :: DIAG_DAVIDSON_QE_ETHR = 1.0d-5
  INTEGER :: IALG_DIAG = 1
  ! 1 => Davidson QE version
  ! 2 => Davidson

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

