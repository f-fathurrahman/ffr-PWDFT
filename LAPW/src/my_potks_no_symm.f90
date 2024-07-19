! A special version of my_potks where potxc is called without symmetrization
! Used for debugging symmetrization routines
!
SUBROUTINE my_potks_no_symm(txc)
  USE m_atoms, ONLY: idxis, natmtot
  USE m_spin, ONLY: spinpol
  USE m_muffin_tins, ONLY: npmt
  USE m_density_pot_xc, ONLY: xcgrad, xctype, wxcmt, wxcir, vxcmt, vxcir, vsir, &
             taumt, tauir, rhomt, rhoir, nosource, msmooth, magmt, magir, &
             exmt, exir, ecmt, ecir, bxcmt, bxcir, vsmt, vclir, vclmt
  USE m_timing, ONLY: timepot
  IMPLICIT NONE 
  ! arguments
  LOGICAL, intent(in) :: txc
  ! local variables
  INTEGER :: is,ias,np
  REAL(8) :: ts0,ts1

  CALL timesec(ts0)

  ! compute the Coulomb potential
  CALL my_potcoul()

  ! meta-GGA variables if required
  IF( (xcgrad==3) .or. (xcgrad==4) ) THEN 
    ! generate the kinetic energy density
    CALL gentau()
    ! compute the Tran-Blaha '09 constant if required
    CALL xc_c_tb09()
  ENDIF 
  
  ! compute the exchange-correlation potential and fields
  !
  ! !!!!! tsh is set to false, symmetrization is not done here !!!!!!!
  !
  IF(txc) CALL my_potxc(.false., xctype, rhomt, rhoir, magmt, magir, taumt, tauir, exmt, &
   exir, ecmt, ecir, vxcmt, vxcir, bxcmt, bxcir, wxcmt, wxcir)

  write(*,*) 'xctype = ', xctype
  !
  write(*,*) 'shape(exmt) = ', shape(exmt)
  write(*,*) 'shape(ecmt) = ', shape(ecir)
  write(*,*) 'shape(vxcmt) = ', shape(vxcmt)
  !
  write(*,*) 'shape(exir) = ', shape(exir)
  write(*,*) 'shape(ecir) = ', shape(ecir)
  write(*,*) 'shape(vxcir) = ', shape(vxcir)
  
  ! optimised effective potential exchange potential
  IF(xctype(1) < 0) CALL oepmain()
  
  ! remove the source term of the exchange-correlation magnetic field if required
  IF(spinpol .and. nosource) CALL projsbf()
  
  ! effective potential from sum of Coulomb and exchange-correlation potentials
  DO ias = 1,natmtot
    is = idxis(ias)
    np = npmt(is)
    vsmt(1:np,ias) = vclmt(1:np,ias) + vxcmt(1:np,ias)
  ENDDO 
  vsir(:) = vclir(:) + vxcir(:)
  
  write(*,*) 'msmooth = ', msmooth
  ! smooth the interstitial potential if required
  CALL rfirsm(msmooth, vsir) ! not needed

  ! generate the effective magnetic fields
  CALL genbs() ! early return for non-spinpol
  
  ! generate the tau-DFT effective potential
  CALL genws() ! only for meta-GGA
  
  CALL timesec(ts1)
  timepot = timepot + ts1 - ts0

  write(*,*) 'sum(vclir) = ', sum(vclir)
  write(*,*) 'sum(vxcir) = ', sum(vxcir)
  write(*,*) 'sum(vsir) = ', sum(vsir)
  write(*,'(1x,A,F18.5)') 'timepot = ', timepot
  
  RETURN 
END SUBROUTINE 