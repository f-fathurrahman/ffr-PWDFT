! !DESCRIPTION:
!   Generates a set of spherical harmonics, $Y_{lm}(\hat{\bf G})$, with angular
!   momenta up to {\tt lmaxo} for the set of ${\bf G}$-vectors.
SUBROUTINE genylmg()
  USE m_gvectors, ONLY: ngvec, vgc, ylmg
  USE m_muffin_tins, ONLY: lmaxo, lmmaxo
  IMPLICIT NONE 
  ! local variables
  INTEGER :: ig
  ! allocate global G-vector spherical harmonic array
  IF(allocated(ylmg)) DEALLOCATE(ylmg)
  ALLOCATE(ylmg(lmmaxo,ngvec))
  DO ig=1,ngvec
    CALL genylmv(lmaxo,vgc(:,ig),ylmg(:,ig))
  ENDDO 
  RETURN 
END SUBROUTINE 