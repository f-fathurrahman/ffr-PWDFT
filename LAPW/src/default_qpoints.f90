SUBROUTINE default_qpoints()
  USE m_qpoints, ONLY: reduceq, ngridq
  IMPLICIT NONE
  reduceq=1
  ngridq(:)=-1
END SUBROUTINE 