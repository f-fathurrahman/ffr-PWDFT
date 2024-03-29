SUBROUTINE rhomagsh()
  ! !USES:
  USE modmain
  ! !DESCRIPTION:
  !   Converts the muffin-tin density and magnetisation from spherical coordinates
  !   to a spherical harmonic expansion. See {\tt rhomagk}.
  IMPLICIT NONE 
  ! local variables
  INTEGER :: idm,is,ias
  INTEGER :: nrc,nrci,npc
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rfmt(:)

  ALLOCATE(rfmt(npcmtmax))
  DO ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
    ! convert the density to spherical harmonics
    CALL dcopy(npc,rhomt(:,ias),1,rfmt,1)
    CALL rfsht(nrc,nrci,rfmt,rhomt(:,ias))
    ! convert magnetisation to spherical harmonics
    DO idm=1,ndmag
      CALL dcopy(npc,magmt(:,ias,idm),1,rfmt,1)
      CALL rfsht(nrc,nrci,rfmt,magmt(:,ias,idm))
    ENDDO 
  ENDDO 
  DEALLOCATE(rfmt)
  RETURN 
END SUBROUTINE 
