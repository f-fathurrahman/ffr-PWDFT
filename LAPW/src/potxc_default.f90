! a simpler version of potxc, all arguments are taken from m_density_pot_xc
! tsh is taken to be .true.
SUBROUTINE potxc_default()
  !
  USE m_atoms, ONLY: natmtot
  USE m_spin, ONLY: spinpol, ncmag
  USE m_muffin_tins, ONLY: npmtmax, npmt, nrmt, nrmti
  USE m_density_pot_xc, ONLY: xctype, wxcmt, wxcir, vxcmt, vxcir, &
             taumt, tauir, rhomt, rhoir, magmt, magir, &
             exmt, exir, ecmt, ecir, bxcmt, bxcir
  !
  IMPLICIT NONE 
  ! arguments
  LOGICAL :: tsh
  ! local variables
  INTEGER :: ias

  tsh = .true.

  ! muffin-tin exchange-correlation potential, field and energy density
  DO ias = 1,natmtot
    CALL potxcmt(tsh, ias, xctype, rhomt, magmt, taumt, exmt, ecmt, vxcmt, bxcmt, wxcmt)
  ENDDO 

  ! interstitial exchange-correlation potential, field and energy density
  CALL potxcir(xctype, rhoir, magir, tauir, exir, ecir, vxcir, bxcir, wxcir)

  ! symmetrise the exchange-correlation potential and magnetic field
  IF(tsh) THEN 
    CALL symrf(nrmt, nrmti, npmt, npmtmax, vxcmt, vxcir)
    IF(spinpol) CALL symrvf(.true., ncmag, nrmt, nrmti, npmt, npmtmax, bxcmt, bxcir)
  ENDIF 
  RETURN 
END SUBROUTINE 
