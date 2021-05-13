SUBROUTINE genrmesh()
  USE m_atoms, ONLY: nspecies, rsp, nrsp, nrspmax, rminsp, rmaxsp
  USE m_muffin_tins, ONLY: nrcmt, rlmt, nrmt, rlcmt, nrmti, nrcmti, rcmt, rmt, &
                     nrcmtmax, nrmtmax, lradstp, lmmaxi, lmmaxo, lmaxo, fracinr, &
                     wprcmt, wrcmt, wrmt, wprmt
! !DESCRIPTION:
!   Generates the coarse and fine radial meshes for each atomic species in the
!   crystal. Also determines which points are in the inner part of the
!   muffin-tin using the value of {\tt fracinr}.

  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,nr,nrc
  INTEGER :: ir,irc,l
  REAL(8) :: t1,t2

  ! estimate the number of radial mesh points to infinity
  nrspmax=1
  DO is=1,nspecies
    ! logarithmic mesh
    t1=log(rmaxsp(is)/rminsp(is))/log(rmt(is)/rminsp(is))
    t2=dble(nrmt(is)-1)*t1
    nrsp(is)=nint(t2)+1
    nrspmax=max(nrspmax,nrsp(is))
  ENDDO 

  ! generate the radial meshes
  IF(allocated(rsp)) DEALLOCATE(rsp)
  ALLOCATE(rsp(nrspmax,nspecies))
  IF(allocated(rlmt)) DEALLOCATE(rlmt)
  ALLOCATE(rlmt(nrmtmax,-lmaxo-1:lmaxo+2,nspecies))
  IF(allocated(wrmt)) DEALLOCATE(wrmt)
  ALLOCATE(wrmt(nrmtmax,nspecies))
  IF(allocated(wprmt)) DEALLOCATE(wprmt)
  ALLOCATE(wprmt(4,nrmtmax,nspecies))
  DO is=1,nspecies
    t1=1.d0/dble(nrmt(is)-1)
    ! logarithmic mesh
    t2=log(rmt(is)/rminsp(is))
    DO ir=1,nrsp(is)
      rsp(ir,is)=rminsp(is)*exp(dble(ir-1)*t1*t2)
    ENDDO 
    ! calculate r^l on the fine radial mesh
    nr=nrmt(is)
    rlmt(1:nr,-1,is)=1.d0/rsp(1:nr,is)
    rlmt(1:nr,0,is)=1.d0
    rlmt(1:nr,1,is)=rsp(1:nr,is)
    DO l=-2,-lmaxo-1,-1
      DO ir=1,nr
        rlmt(ir,l,is)=rlmt(ir,l+1,is)/rsp(ir,is)
      ENDDO 
    ENDDO 
    DO l=2,lmaxo+2
      DO ir=1,nr
        rlmt(ir,l,is)=rlmt(ir,l-1,is)*rsp(ir,is)
      ENDDO 
    ENDDO 
    ! determine the weights for spline integration on the fine radial mesh
    CALL wsplint(nr,rsp(:,is),wrmt(:,is))
    ! multiply by r^2
    wrmt(1:nr,is)=wrmt(1:nr,is)*rlmt(1:nr,2,is)
    ! determine the weights for partial integration on fine radial mesh
    CALL wsplintp(nr,rsp(:,is),wprmt(:,:,is))
  ENDDO 

! determine the fraction of the muffin-tin radius which defines the inner part
IF(fracinr.lt.0.d0) fracinr=sqrt(dble(lmmaxi)/dble(lmmaxo))

! set up the coarse radial meshes and find the inner part of the muffin-tin
! where rho is calculated with lmaxi
IF(allocated(rcmt)) DEALLOCATE(rcmt)
ALLOCATE(rcmt(nrcmtmax,nspecies))
IF(allocated(rlcmt)) DEALLOCATE(rlcmt)
ALLOCATE(rlcmt(nrcmtmax,-lmaxo-1:lmaxo+2,nspecies))
IF(allocated(wrcmt)) DEALLOCATE(wrcmt)
ALLOCATE(wrcmt(nrcmtmax,nspecies))
IF(allocated(wprcmt)) DEALLOCATE(wprcmt)
ALLOCATE(wprcmt(4,nrcmtmax,nspecies))
DO is=1,nspecies
  t1=fracinr*rmt(is)
  nrmti(is)=1
  nrcmti(is)=1
  irc=0
  DO ir=1,nrmt(is),lradstp
    irc=irc+1
    rcmt(irc,is)=rsp(ir,is)
    IF(rsp(ir,is).lt.t1) THEN 
      nrmti(is)=ir
      nrcmti(is)=irc
    ENDIF 
  ENDDO 
! store r^l on the coarse radial mesh
  DO l=-lmaxo-1,lmaxo+2
    irc=0
    DO ir=1,nrmt(is),lradstp
      irc=irc+1
      rlcmt(irc,l,is)=rlmt(ir,l,is)
    ENDDO 
  ENDDO 
! determine the weights for spline integration on the coarse radial mesh
  nrc=nrcmt(is)
  CALL wsplint(nrc,rcmt(:,is),wrcmt(:,is))
! multiply by r^2
  wrcmt(1:nrc,is)=wrcmt(1:nrc,is)*rlcmt(1:nrc,2,is)
! determine the weights for partial integration on coarse radial mesh
  CALL wsplintp(nrc,rcmt(:,is),wprcmt(:,:,is))
ENDDO 

RETURN 
END SUBROUTINE 
!EOC

