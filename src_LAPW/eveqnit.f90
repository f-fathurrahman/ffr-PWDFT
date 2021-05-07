SUBROUTINE eveqnit(nmatp,ngp,igpig,vpl,vgpl,vgpc,apwalm,evalfv,evecfv)
use modmain
IMPLICIT NONE 
! arguments
INTEGER, intent(in) :: nmatp,ngp
INTEGER, intent(in) :: igpig(ngkmax)
REAL(8), intent(in) :: vpl(3)
REAL(8), intent(in) :: vgpl(3,ngkmax),vgpc(3,ngkmax)
COMPLEX(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
REAL(8), intent(out) :: evalfv(nstfv)
COMPLEX(8), intent(out) :: evecfv(nmatmax,nstfv)
  ! local variables
  INTEGER :: n2,ns,ist,ias
  INTEGER :: it,i
  INTEGER :: lwork,info
  REAL(8) :: rmax,t1
  REAL(8) :: ts1,ts0
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: w(:),rwork(:)
  COMPLEX(8), ALLOCATABLE :: h(:,:),o(:,:),hv(:,:),ov(:,:)
  COMPLEX(8), ALLOCATABLE :: u(:,:),hu(:,:),ou(:,:)
  COMPLEX(8), ALLOCATABLE :: hs(:,:),os(:,:),work(:)
  ! external functions
  COMPLEX(8) zdotc
  external zdotc
  n2=2*nmatp
  ns=2*nstfv

  IF(iscl.ge.2) THEN 
  ! read in the eigenvalues/vectors from file
    CALL getevalfv(filext,0,vpl,evalfv)
    CALL getevecfv(filext,0,vpl,vgpl,evecfv)
  ELSE 
  ! initialise the eigenvectors to canonical basis vectors
    evecfv(1:nmatp,:)=0.d0
    DO ist=1,nstfv
      evecfv(ist,ist)=1.d0
    ENDDO 
  ENDIF 
  ! compute Hamiltonian and overlap matrices
  CALL timesec(ts0)
  ALLOCATE(h(nmatp,nmatp),o(nmatp,nmatp))
  ! Hamiltonian
  DO i=1,nmatp
    h(1:i,i)=0.d0
  ENDDO 
  DO ias=1,natmtot
    CALL hmlaa(tefvr,ias,ngp,apwalm(:,:,:,ias),nmatp,h)
    CALL hmlalo(ias,ngp,apwalm(:,:,:,ias),nmatp,h)
    CALL hmllolo(ias,ngp,nmatp,h)
  ENDDO 
  CALL hmlistl(ngp,igpig,vgpc,nmatp,h)
  ! overlap
  DO i=1,nmatp
    o(1:i,i)=0.d0
  ENDDO 
  DO ias=1,natmtot
    CALL olpaa(tefvr,ias,ngp,apwalm(:,:,:,ias),nmatp,o)
    CALL olpalo(ias,ngp,apwalm(:,:,:,ias),nmatp,o)
    CALL olplolo(ias,ngp,nmatp,o)
  ENDDO 
  CALL olpistl(ngp,igpig,nmatp,o)

  CALL timesec(ts1)

  timemat=timemat+ts1-ts0
  CALL timesec(ts0)
  ALLOCATE(w(ns),rwork(3*ns))
  ALLOCATE(hv(nmatp,nstfv),ov(nmatp,nstfv))
  ALLOCATE(u(nmatp,nstfv),hu(nmatp,nstfv),ou(nmatp,nstfv))
  ALLOCATE(hs(ns,ns),os(ns,ns))
  lwork=2*ns
  ALLOCATE(work(lwork))

  ! start iteration loop
  DO it=1,maxitefv
    rmax=0.d0
    DO ist=1,nstfv
      ! operate with O on the current eigenvector
      CALL zhemv('U',nmatp,zone,o,nmatp,evecfv(:,ist),1,zzero,ov(:,ist),1)
      ! normalise the eigenvector
      t1=dble(zdotc(nmatp,evecfv(:,ist),1,ov(:,ist),1))
      IF(t1.gt.0.d0) THEN 
        t1=1.d0/sqrt(t1)
        CALL dscal(n2,t1,evecfv(:,ist),1)
        CALL dscal(n2,t1,ov(:,ist),1)
      ENDIF 
      ! operate with H on the current eigenvector
      CALL zhemv('U',nmatp,zone,h,nmatp,evecfv(:,ist),1,zzero,hv(:,ist),1)
      ! estimate the eigenvalue
      t1=dble(zdotc(nmatp,evecfv(:,ist),1,hv(:,ist),1))
      IF((iscl.le.1).and.(it.eq.1)) THEN 
        evalfv(ist)=t1
      ELSE 
        evalfv(ist)=(1.d0-befvit)*evalfv(ist)+befvit*t1
      ENDIF 
      ! compute the residual |u> = (H - eO)|v>
      CALL zcopy(nmatp,hv(:,ist),1,u(:,ist),1)
      CALL daxpy(n2,-evalfv(ist),ov(:,ist),1,u(:,ist),1)
      ! apply the overlap matrix to the residual
      CALL zhemv('U',nmatp,zone,o,nmatp,u(:,ist),1,zzero,ou(:,ist),1)
      ! compute the overlap of the residual with itself
      t1=dble(zdotc(nmatp,u(:,ist),1,ou(:,ist),1))
      IF(t1.gt.rmax) rmax=t1
       ! normalise the residual
      IF(t1.gt.0.d0) THEN 
        t1=1.d0/sqrt(t1)
        CALL dscal(n2,t1,u(:,ist),1)
        CALL dscal(n2,t1,ou(:,ist),1)
      ENDIF 
      ! apply the Hamiltonian matrix to the residual
      CALL zhemv('U',nmatp,zone,h,nmatp,u(:,ist),1,zzero,hu(:,ist),1)
    ENDDO 
  ! compute the Hamiltonian and overlap matrices in the subspace formed by the
  ! eigenvectors and their residuals
    DO ist=1,nstfv
      CALL zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,hv(:,ist),1,zzero, &
       hs(1,ist),1)
      CALL zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,hu(:,ist),1,zzero, &
       hs(1,nstfv+ist),1)
      CALL zgemv('C',nmatp,nstfv,zone,u,nmatp,hu(:,ist),1,zzero, &
       hs(nstfv+1,nstfv+ist),1)
    ENDDO 
  
    DO ist=1,nstfv
      CALL zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,ov(:,ist),1,zzero, &
       os(1,ist),1)
      CALL zgemv('C',nmatp,nstfv,zone,evecfv,nmatmax,ou(:,ist),1,zzero, &
       os(1,nstfv+ist),1)
      CALL zgemv('C',nmatp,nstfv,zone,u,nmatp,ou(:,ist),1,zzero, &
       os(nstfv+1,nstfv+ist),1)
    ENDDO 
  ! solve the generalised eigenvalue problem in the subspace
    CALL zhegv(1,'V','U',ns,hs,ns,os,ns,w,work,lwork,rwork,info)
    IF(info.ne.0) exit
    ! construct the new eigenvectors
    DO ist=1,nstfv
      CALL zgemv('N',nmatp,nstfv,zone,evecfv,nmatmax,hs(1,ist),1,zzero, &
       ov(:,ist),1)
      CALL zgemv('N',nmatp,nstfv,zone,u,nmatp,hs(nstfv+1,ist),1,zone,ov(:,ist),1)
    ENDDO 
  
    DO ist=1,nstfv
      CALL zcopy(nmatp,ov(:,ist),1,evecfv(:,ist),1)
    ENDDO 
    
    ! check for convergence
    rmax=sqrt(abs(rmax)/dble(nmatp))
    IF((it.ge.minitefv).and.(rmax.lt.epsefvit)) exit
    ! end iteration loop
  ENDDO 
  DEALLOCATE(w,rwork,h,o,hv,ov)
  DEALLOCATE(u,hu,ou,hs,os,work)
  CALL timesec(ts1)
  timefv=timefv+ts1-ts0
  RETURN 
END SUBROUTINE 

