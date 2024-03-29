SUBROUTINE curden(afield)
  use modmain
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(in) :: afield(3)
  ! local variables
  INTEGER ik,is,ias,ispn
  INTEGER nr,nri,iro,np
  INTEGER nrc,nrci,npc
  INTEGER ir,n,i,nthd
  REAL(8) ca,t1
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rfmt(:)
  ! external functions
  REAL(8) rfint
  external rfint

  ! coupling constant of the external A-field (1/c)
  ca=1.d0/solsc
  ! set the current density to zero
  DO i=1,3
    DO ias=1,natmtot
      is=idxis(ias)
      cdmt(1:npcmt(is),ias,i)=0.d0
    ENDDO 
  ENDDO 
  cdir(:,:)=0.d0

  ! current density cannot be computed if wavefunctions DO not exist
  IF(iscl <= 0) RETURN 
  DO ik=1,nkpt
    CALL curdenk(ik)
  ENDDO 

  ! convert muffin-tin current density to spherical harmonics
  ALLOCATE(rfmt(npcmtmax))
  DO ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
    DO i=1,3
      rfmt(1:npc)=cdmt(1:npc,ias,i)
      CALL rfsht(nrc,nrci,rfmt,cdmt(:,ias,i))
    ENDDO 
  ENDDO 
  DEALLOCATE(rfmt)

  ! symmetrise the current density
  CALL symrvf(.false.,.true.,nrcmt,nrcmti,npcmt,npmtmax,cdmt,cdir)

  ! convert the current density from a coarse to a fine radial mesh
  DO i=1,3
    CALL rfmtctof(cdmt(:,:,i))
  ENDDO 

  ! add vector potential contribution to make current gauge invariant
  ALLOCATE(rfmt(npmtmax))
  DO ias=1,natmtot
    is=idxis(ias)
    nr=nrmt(is)
    nri=nrmti(is)
    np=npmt(is)
    iro=nri+1
    ! remove the core density from the muffin-tin density
    CALL dcopy(np,rhomt(:,ias),1,rfmt,1)
    DO ispn=1,nspncr
      i=1
      DO ir=1,nri
        rfmt(i)=rfmt(i)-rhocr(ir,ias,ispn)
        i=i+lmmaxi
      ENDDO 
      DO ir=iro,nr
        rfmt(i)=rfmt(i)-rhocr(ir,ias,ispn)
        i=i+lmmaxo
      ENDDO 
    ENDDO 
    DO i=1,3
      t1=-ca*afield(i)
      CALL daxpy(np,t1,rfmt,1,cdmt(:,ias,i),1)
    ENDDO 
  ENDDO 
  DEALLOCATE(rfmt)
  DO i=1,3
    t1=-ca*afield(i)
    CALL daxpy(ngtot,t1,rhoir,1,cdir(:,i),1)
  ENDDO 
  ! compute the total current in the unit cell
  DO i=1,3
    curtot(i)=rfint(cdmt(:,:,i),cdir(:,i))
  ENDDO 
  ! total current magnitude
  curtotm=sqrt(curtot(1)**2+curtot(2)**2+curtot(3)**2)
  RETURN 
END SUBROUTINE 

