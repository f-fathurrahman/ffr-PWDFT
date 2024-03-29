SUBROUTINE genbs()
  use modmain
  IMPLICIT NONE 
  ! local variables
  INTEGER :: idm,is,ia,ias
  INTEGER :: nrc,nrci,npc
  REAL(8) :: cb,t1
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rfmt(:)

  IF(.not. spinpol) RETURN 
  
  ! coupling constant of the external field (g_e/4c)
  cb = gfacte/(4.d0*solsc)

  !------------------------------------!
  !     muffin-tin Kohn-Sham field     !
  !------------------------------------!
  ALLOCATE(rfmt(npcmtmax))
  DO ias=1,natmtot
    is=idxis(ias)
    ia=idxia(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
    ! exchange-correlation magnetic field in spherical coordinates
    DO idm=1,ndmag
      CALL rfmtftoc(nrc,nrci,bxcmt(:,ias,idm),rfmt)
      CALL rbsht(nrc,nrci,rfmt,bsmt(:,ias,idm))
    ENDDO 
    ! add the external magnetic field
    t1=cb*(bfcmt(3,ia,is)+bfieldc(3))
    bsmt(1:npc,ias,ndmag)=bsmt(1:npc,ias,ndmag)+t1
    IF(ncmag) THEN 
      DO idm=1,2
        t1=cb*(bfcmt(idm,ia,is)+bfieldc(idm))
        bsmt(1:npc,ias,idm)=bsmt(1:npc,ias,idm)+t1
      ENDDO 
    ENDIF 
  ENDDO 
  DEALLOCATE(rfmt)

  !-----------------------------------------------!
  !     interstitial Kohn-Sham magnetic field     !
  !-----------------------------------------------!
  DO idm=1,ndmag
    IF(ncmag) THEN 
      t1 = cb*bfieldc(idm)
    ELSE 
      t1 = cb*bfieldc(3)
    ENDIF 
    ! multiply by characteristic function
    bsir(:,idm)=(bxcir(:,idm)+t1)*cfunir(:)
  ENDDO 
  ! add the magnetic dipole field if required
  IF(tbdip) CALL bdipole()
  RETURN 
END SUBROUTINE 

