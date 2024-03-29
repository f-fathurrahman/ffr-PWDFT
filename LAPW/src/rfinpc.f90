REAL(8) FUNCTION rfinpc(ld,rfmt1,rfir1,rfmt2,rfir2)
  use modmain
  IMPLICIT NONE 
  INTEGER, intent(in) :: ld
  REAL(8), intent(in) :: rfmt1(ld,natmtot),rfir1(ngtot)
  REAL(8), intent(in) :: rfmt2(ld,natmtot),rfir2(ngtot)
  ! local variables
  INTEGER is,ias,ir
  ! external functions
  REAL(8) rfmtinp
  external rfmtinp
  ! interstitial contribution
  rfinpc=0.d0
  DO ir=1,ngtot
    rfinpc=rfinpc+rfir1(ir)*rfir2(ir)*cfunir(ir)
  ENDDO 
  rfinpc=rfinpc*omega/dble(ngtot)
  ! muffin-tin contribution
  DO ias=1,natmtot
    is=idxis(ias)
    rfinpc=rfinpc+rfmtinp(nrcmt(is),nrcmti(is),wrcmt(:,is),rfmt1(:,ias), &
     rfmt2(:,ias))
  ENDDO 

  RETURN 
END FUNCTION 