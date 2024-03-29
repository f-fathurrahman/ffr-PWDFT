SUBROUTINE gentau()
  USE modmain
  IMPLICIT NONE 
  ! local variables
  INTEGER ik,ispn,is,ias
  INTEGER np,npc,n,nthd
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rfmt(:,:),rfir(:)
  REAL(8), ALLOCATABLE :: rvfmt(:,:,:),rvfir(:,:)

  ! set the kinetic energy density to zero
  taumt(:,:,:)=0.d0
  tauir(:,:)=0.d0

  ! tau cannot be computed if wavefunctions DO not exist
  IF(iscl <= 0) RETURN 
  DO ik=1,nkpt
    CALL gentauk(ik)
  ENDDO 
  ALLOCATE(rfmt(npmtmax,natmtot))
  ! convert taumt to spherical harmonics
  DO ispn=1,nspinor
    DO ias=1,natmtot
      is=idxis(ias)
      CALL dcopy(npcmt(is),taumt(:,ias,ispn),1,rfmt,1)
      CALL rfsht(nrcmt(is),nrcmti(is),rfmt,taumt(:,ias,ispn))
    ENDDO 
  ENDDO 

  ! symmetrise tau
  IF(spinpol) THEN 
    ! spin-polarised case: convert to scalar-vector form
    ALLOCATE(rfir(ngtot))
    ALLOCATE(rvfmt(npmtmax,natmtot,ndmag))
    ALLOCATE(rvfir(ngtot,ndmag))
    DO ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      rfmt(1:npc,ias)=taumt(1:npc,ias,1)+taumt(1:npc,ias,2)
      rvfmt(1:npc,ias,1:ndmag-1)=0.d0
      rvfmt(1:npc,ias,ndmag)=taumt(1:npc,ias,1)-taumt(1:npc,ias,2)
    ENDDO 
    rfir(:)=tauir(:,1)+tauir(:,2)
    rvfir(:,1:ndmag-1)=0.d0
    rvfir(:,ndmag)=tauir(:,1)-tauir(:,2)
    CALL symrf(nrcmt,nrcmti,npcmt,npmtmax,rfmt,rfir)
    CALL symrvf(.true.,ncmag,nrcmt,nrcmti,npcmt,npmtmax,rvfmt,rvfir)
    DO ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      taumt(1:npc,ias,1)=0.5d0*(rfmt(1:npc,ias)+rvfmt(1:npc,ias,ndmag))
      taumt(1:npc,ias,2)=0.5d0*(rfmt(1:npc,ias)-rvfmt(1:npc,ias,ndmag))
    ENDDO 
    tauir(:,1)=0.5d0*(rfir(:)+rvfir(:,ndmag))
    tauir(:,2)=0.5d0*(rfir(:)-rvfir(:,ndmag))
    DEALLOCATE(rfir,rvfmt,rvfir)
  ELSE 
    ! spin-unpolarised case
    CALL symrf(nrcmt,nrcmti,npcmt,npmtmax,taumt,tauir)
  ENDIF 

  ! convert taumt from a coarse to a fine radial mesh
  DO ispn=1,nspinor
    CALL rfmtctof(taumt(:,:,ispn))
  ENDDO 

  ! generate the core kinetic energy density
  CALL gentaucr()
  DO ispn=1,nspinor
    DO ias=1,natmtot
      is=idxis(ias)
      np=npmt(is)
      ! add the core contribution
      taumt(1:np,ias,ispn)=taumt(1:np,ias,ispn)+taucr(1:np,ias,ispn)
      ! zero tau on the inner part of the muffin-tin
      taumt(1:npmti(is),ias,ispn)=0.d0
    ENDDO 
  ENDDO 
  DEALLOCATE(rfmt)
  RETURN 
END SUBROUTINE 

