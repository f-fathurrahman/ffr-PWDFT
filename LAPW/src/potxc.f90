SUBROUTINE potxc(tsh,xctype_,rhomt_,rhoir_,magmt_,magir_,taumt_,tauir_,exmt_, &
 exir_,ecmt_,ecir_,vxcmt_,vxcir_,bxcmt_,bxcir_,wxcmt_,wxcir_)
  !
  USE m_atoms, ONLY: natmtot
  USE m_gvectors, ONLY: ngtot
  USE m_spin, ONLY: ndmag, nspinor, spinpol, ncmag
  USE m_muffin_tins, ONLY: npmtmax, npmt, nrmt, nrmti
  IMPLICIT NONE 
  ! arguments
  LOGICAL, intent(in) :: tsh
  INTEGER, intent(in) :: xctype_(3)
  REAL(8), intent(in) :: rhomt_(npmtmax,natmtot),rhoir_(ngtot)
  REAL(8), intent(in) :: magmt_(npmtmax,natmtot,ndmag),magir_(ngtot,ndmag)
  REAL(8), intent(in) :: taumt_(npmtmax,natmtot,nspinor),tauir_(ngtot,nspinor)
  REAL(8), intent(out) :: exmt_(npmtmax,natmtot),exir_(ngtot)
  REAL(8), intent(out) :: ecmt_(npmtmax,natmtot),ecir_(ngtot)
  REAL(8), intent(out) :: vxcmt_(npmtmax,natmtot),vxcir_(ngtot)
  REAL(8), intent(out) :: bxcmt_(npmtmax,natmtot,ndmag),bxcir_(ngtot,ndmag)
  REAL(8), intent(out) :: wxcmt_(npmtmax,natmtot),wxcir_(ngtot)
  ! local variables
  INTEGER :: ias

  ! muffin-tin exchange-correlation potential, field and energy density
  DO ias=1,natmtot
    CALL potxcmt(tsh,ias,xctype_,rhomt_,magmt_,taumt_,exmt_,ecmt_,vxcmt_,bxcmt_, &
     wxcmt_)
  ENDDO 

  ! interstitial exchange-correlation potential, field and energy density
  CALL potxcir(xctype_,rhoir_,magir_,tauir_,exir_,ecir_,vxcir_,bxcir_,wxcir_)

  ! symmetrise the exchange-correlation potential and magnetic field
  IF(tsh) THEN 
    CALL symrf(nrmt,nrmti,npmt,npmtmax,vxcmt_,vxcir_)
    IF(spinpol) CALL symrvf(.true.,ncmag,nrmt,nrmti,npmt,npmtmax,bxcmt_,bxcir_)
  ENDIF 
  RETURN 
END SUBROUTINE 
