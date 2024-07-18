SUBROUTINE my_potxc(tsh, xctype_, &
  rhomt_, rhoir_, magmt_, magir_, taumt_,tauir_, &
  exmt_, exir_, ecmt_, ecir_, &
  vxcmt_, vxcir_, bxcmt_, bxcir_, wxcmt_, wxcir_ )
  !
  USE m_atoms, ONLY: natmtot
  USE m_gvectors, ONLY: ngtot
  USE m_spin, ONLY: ndmag, nspinor, spinpol, ncmag
  USE m_muffin_tins, ONLY: npmtmax, npmt, nrmt, nrmti
  IMPLICIT NONE 
  ! arguments
  LOGICAL, intent(in) :: tsh   ! FIXME: always true?
  INTEGER, intent(in) :: xctype_(3)
  REAL(8), intent(in) :: rhomt_(npmtmax,natmtot), rhoir_(ngtot)
  REAL(8), intent(in) :: magmt_(npmtmax,natmtot,ndmag), magir_(ngtot,ndmag)
  REAL(8), intent(in) :: taumt_(npmtmax,natmtot,nspinor), tauir_(ngtot,nspinor)
  !
  REAL(8), intent(out) :: exmt_(npmtmax,natmtot), exir_(ngtot)
  REAL(8), intent(out) :: ecmt_(npmtmax,natmtot), ecir_(ngtot)
  REAL(8), intent(out) :: vxcmt_(npmtmax,natmtot), vxcir_(ngtot)
  REAL(8), intent(out) :: bxcmt_(npmtmax,natmtot,ndmag), bxcir_(ngtot,ndmag)
  REAL(8), intent(out) :: wxcmt_(npmtmax,natmtot), wxcir_(ngtot)
  ! local variables
  INTEGER :: ias

  write(*,*)
  write(*,*) 'my_potxc: shape(rhomt) = ', shape(rhomt_)
  write(*,*) 'my_potxc: sum(rhomt) = ', sum(rhomt_)
  !
  write(*,*) 'my_potxc: shape(magmt) = ', shape(magmt_)
  write(*,*) 'my_potxc: sum(magmt) = ', sum(magmt_)
  !
  !write(*,*) 'my_potxc: shape(taumt) = ', shape(taumt_)
  ! XXX: Will result in segfault if not a meta-GGA calc
  !write(*,*) 'my_potxc: sum(taumt) = ', sum(taumt_)

  ! muffin-tin exchange-correlation potential, field and energy density
  ! This is called for each atom
  DO ias = 1,natmtot
    !CALL my_potxcmt( tsh, ias, xctype_, &
    !  rhomt_, magmt_, taumt_, exmt_, ecmt_, vxcmt_, bxcmt_, wxcmt_ )
    CALL my_potxcmt_lda_nospin( tsh, ias, xctype_, rhomt_, exmt_, ecmt_, vxcmt_ )
  ENDDO 

  ! interstitial exchange-correlation potential, field and energy density
  !CALL my_potxcir(xctype_,rhoir_,magir_,tauir_,exir_,ecir_,vxcir_,bxcir_,wxcir_)
  CALL my_potxcir_lda_nospin(xctype_, rhoir_, exir_, ecir_, vxcir_)
  
  write(*,*) 'sum exmt_ = ', sum(exmt_)
  write(*,*) 'sum ecmt_ = ', sum(ecmt_)
  write(*,*) 'sum exir_ = ', sum(exir_)
  write(*,*) 'sum ecir_ = ', sum(ecir_)

  ! symmetrize the exchange-correlation potential and magnetic field
  write(*,*) 'Before symmetrizing: sum vxcmt = ', sum(vxcmt_)
  write(*,*) 'Before symmetrizing: sum vxcir = ', sum(vxcir_)
  IF(tsh) THEN 
    write(*,*) 'my_potxc: symmetrizing Vxc'
    CALL my_symrf(nrmt, nrmti, npmt, npmtmax, vxcmt_, vxcir_)
    IF(spinpol) CALL symrvf(.true., ncmag, nrmt, nrmti, npmt, npmtmax, bxcmt_, bxcir_)
  ENDIF 
  write(*,*) 'After symmetrizing: sum vxcmt = ', sum(vxcmt_)
  write(*,*) 'After symmetrizing: sum vxcir = ', sum(vxcir_)


  stop 'ffr 45 in my_potxc'

  RETURN 
END SUBROUTINE 
