SUBROUTINE genws()
  use modmain
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ias
  INTEGER :: nrc,nrci
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rfmt(:)
  
  IF(xcgrad /= 4) RETURN 
  
  ! muffin-tin effective tau-DFT potential
  ALLOCATE(rfmt(npcmtmax))
  DO ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    ! convert to coarse radial mesh and spherical coordinates
    CALL rfmtftoc(nrc,nrci,wxcmt(:,ias),rfmt)
    CALL rbsht(nrc,nrci,rfmt,wsmt(:,ias))
  ENDDO 
  DEALLOCATE(rfmt)

  ! interstitial tau-DFT potential
  wsir(:) = wxcir(:)*cfunir(:)
  RETURN 
END SUBROUTINE 