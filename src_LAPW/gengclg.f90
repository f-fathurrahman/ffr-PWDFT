SUBROUTINE gengclg()
  USE m_gvectors, ONLY: gclg, ngvec, gc
  USE m_constants, ONLY: fourpi
  IMPLICIT NONE 
  ! local variables
  IF(allocated(gclg)) DEALLOCATE(gclg)
  ALLOCATE(gclg(ngvec))
  gclg(1)=0.d0
  gclg(2:ngvec)=fourpi/gc(2:ngvec)**2
  RETURN 
END SUBROUTINE 

