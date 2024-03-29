SUBROUTINE gensocfr()
  USE modmain
  IMPLICIT NONE 
  ! local variables
  INTEGER is,ias,nthd
  INTEGER nr,nri,ir,irc
  REAL(8) cso,rm
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: vr(:),dvr(:)
  IF(.not.spinorb) RETURN 
  ! coefficient of spin-orbit coupling
  cso=socscf/(4.d0*solsc**2)
  ALLOCATE(vr(nrmtmax),dvr(nrmtmax))
  DO ias=1,natmtot
    is=idxis(ias)
    nr=nrmt(is)
    nri=nrmti(is)
  ! radial derivative of the spherical part of the Kohn-Sham potential
    CALL rfmtlm(1,nr,nri,vsmt(:,ias),vr)
    vr(1:nr)=vr(1:nr)*y00
    CALL fderiv(1,nr,rlmt(:,1,is),vr,dvr)
    irc=0
    DO ir=1,nr,lradstp
      irc=irc+1
      rm=1.d0-2.d0*cso*vr(ir)
      socfr(irc,ias)=cso*dvr(ir)/(rsp(ir,is)*rm**2)
    ENDDO 
  ENDDO 
  DEALLOCATE(vr,dvr)
  RETURN 
END SUBROUTINE 