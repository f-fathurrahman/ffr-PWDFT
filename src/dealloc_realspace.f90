SUBROUTINE dealloc_realspace()

  USE m_realspace
  IMPLICIT NONE 

  IF( allocated(rgrid) ) DEALLOCATE( rgrid )

END SUBROUTINE 

