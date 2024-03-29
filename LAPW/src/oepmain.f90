SUBROUTINE oepmain()
  USE modmain
  IMPLICIT NONE 
  ! local variables
  INTEGER :: ik,idm,is,ias
  INTEGER :: nrc,nrci,np,npc
  INTEGER :: n,it
  REAL(8) :: t1
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: dvxmt(:,:),dvxir(:)
  REAL(8), ALLOCATABLE :: dbxmt(:,:,:),dbxir(:,:)
  REAL(8), ALLOCATABLE :: rfmt1(:,:),rfmt2(:),rfir(:)
  REAL(8), ALLOCATABLE :: rvfmt(:,:,:),rvfir(:,:)
  COMPLEX(8), ALLOCATABLE :: vclcv(:,:,:,:),vclvv(:,:,:)
  ! external functions
  REAL(8) rfinpc
  external rfinpc

  ! initialise the OEP exchange potential
  IF(iscl <= 0) THEN 
    CALL initoep()
    RETURN 
  ENDIF 

  ! calculate Coulomb matrix elements
  ALLOCATE(vclcv(ncrmax,natmtot,nstsv,nkpt),vclvv(nstsv,nstsv,nkpt))
  CALL oepvcl(vclcv,vclvv)
  
  ! allocate local arrays
  ALLOCATE(dvxmt(npcmtmax,natmtot),dvxir(ngtot))
  ALLOCATE(rfmt1(npmtmax,natmtot),rfir(ngtot))
  IF(spinpol) THEN 
    ALLOCATE(dbxmt(npcmtmax,natmtot,ndmag),dbxir(ngtot,ndmag))
    ALLOCATE(rvfmt(npmtmax,natmtot,ndmag),rvfir(ngtot,ndmag))
  ENDIF 

  !------------------------------!
  !     start iteration loop     !
  !------------------------------!
  DO it=1,maxitoep
  
    IF( mod(it,10)==0 ) THEN 
      WRITE(*,'("Info(oepmain): done ",I4," iterations of ",I4)') it, maxitoep
    ENDIF 
    ! zero the residuals
    dvxmt(:,:)=0.d0
    dvxir(:)=0.d0
    IF(spinpol) THEN 
      dbxmt(:,:,:)=0.d0
      dbxir(:,:)=0.d0
    ENDIF 
  
    ! calculate the k-dependent residuals
    DO ik=1,nkpt
      CALL oepresk(ik,vclcv,vclvv,dvxmt,dvxir,dbxmt,dbxir)
    ENDDO 
  
    ! convert muffin-tin residuals to spherical harmonics
    DO ias=1,natmtot
      is=idxis(ias)
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      CALL rfsht(nrc,nrci,dvxmt(:,ias),rfmt1(:,ias))
      DO idm=1,ndmag
        CALL rfsht(nrc,nrci,dbxmt(:,ias,idm),rvfmt(:,ias,idm))
      ENDDO 
    ENDDO 
  
    ! symmetrise the residuals
    CALL symrf(nrcmt,nrcmti,npcmt,npmtmax,rfmt1,dvxir)
    IF(spinpol) CALL symrvf(.true.,ncmag,nrcmt,nrcmti,npcmt,npmtmax,rvfmt,dbxir)
    
    ! magnitude of residuals
    resoep=sqrt(abs(rfinpc(npmtmax,rfmt1,dvxir,rfmt1,dvxir)))
    DO idm=1,ndmag
      t1=rfinpc(npmtmax,rvfmt(:,:,idm),dbxir(:,idm),rvfmt(:,:,idm),dbxir(:,idm))
      resoep=resoep+sqrt(abs(t1))
    ENDDO 
    resoep=resoep/omega
    
    ! update exchange potential and magnetic field
    ALLOCATE(rfmt2(npcmtmax))
    DO ias=1,natmtot
      is=idxis(ias)
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
      ! convert residual to spherical coordinates
      CALL rbsht(nrc,nrci,rfmt1(:,ias),rfmt2)
      ! subtract from exchange potential
      vxmt(1:npc,ias)=vxmt(1:npc,ias)-tauoep*rfmt2(1:npc)
      ! repeat for exchange magnetic field
      DO idm=1,ndmag
        CALL rbsht(nrc,nrci,rvfmt(:,ias,idm),rfmt2)
        bxmt(1:npc,ias,idm)=bxmt(1:npc,ias,idm)-tauoep*rfmt2(1:npc)
      ENDDO 
    ENDDO 
    DEALLOCATE(rfmt2)
  
    vxir(:)=vxir(:)-tauoep*dvxir(:)
    DO idm=1,ndmag
      bxir(:,idm)=bxir(:,idm)-tauoep*dbxir(:,idm)
    ENDDO 

  ENDDO ! end iteration loop

  ! convert the exchange potential and field to spherical harmonics
  DO ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    CALL rfsht(nrc,nrci,vxmt(:,ias),rfmt1(:,ias))
    DO idm=1,ndmag
      CALL rfsht(nrc,nrci,bxmt(:,ias,idm),rvfmt(:,ias,idm))
    ENDDO 
  ENDDO 
  
  ! convert potential and field from a coarse to a fine radial mesh
  CALL rfmtctof(rfmt1)
  DO idm=1,ndmag
    CALL rfmtctof(rvfmt(:,:,idm))
  ENDDO 
  
  ! add to existing (density derived) correlation potential and field
  DO ias=1,natmtot
    is=idxis(ias)
    np=npmt(is)
    vxcmt(1:np,ias)=vxcmt(1:np,ias)+rfmt1(1:np,ias)
    DO idm=1,ndmag
      bxcmt(1:np,ias,idm)=bxcmt(1:np,ias,idm)+rvfmt(1:np,ias,idm)
    ENDDO 
  ENDDO 
  vxcir(:)=vxcir(:)+vxir(:)
  DO idm=1,ndmag
    bxcir(:,idm)=bxcir(:,idm)+bxir(:,idm)
  ENDDO 
  ! symmetrise the exchange potential and field
  CALL symrf(nrmt,nrmti,npmt,npmtmax,vxcmt,vxcir)
  IF(spinpol) CALL symrvf(.true.,ncmag,nrmt,nrmti,npmt,npmtmax,bxcmt,bxcir)
  DEALLOCATE(rfmt1,rfir,vclcv,vclvv)
  DEALLOCATE(dvxmt,dvxir)
  IF(spinpol) THEN 
    DEALLOCATE(rvfmt,rvfir)
    DEALLOCATE(dbxmt,dbxir)
  ENDIF 
  ! set the constant part of the exchange potential equal to that of LDA/GGA
  CALL rfint0(vxc0,vxcmt,vxcir)
  RETURN 
END SUBROUTINE 