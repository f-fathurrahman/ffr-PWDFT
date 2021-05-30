PROGRAM main
  IMPLICIT NONE 

  ! read input and initialize some variables
  CALL read_input()
  CALL init0()
  CALL init1()
  
  CALL info_crystal()
  CALL info_symmetry()
  CALL writesym()
  CALL info_gvectors()
  CALL info_muffin_tins()

  CALL rhoinit()
  CALL debug_potks(.true.)

END PROGRAM


SUBROUTINE debug_potks(txc)
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
  IF((xcgrad.eq.3).or.(xcgrad.eq.4)) THEN 
    ! generate the kinetic energy density
    CALL gentau()
    ! compute the Tran-Blaha '09 constant if required
    CALL xc_c_tb09()
  ENDIF 
  
  ! compute the exchange-correlation potential and fields
  IF(txc) CALL my_potxc(.true.,xctype,rhomt,rhoir,magmt,magir,taumt,tauir,exmt, &
   exir,ecmt,ecir,vxcmt,vxcir,bxcmt,bxcir,wxcmt,wxcir)

  write(*,*) 'xctype = ', xctype
  write(*,*) 'shape(exmt) = ', shape(exmt)
  write(*,*) 'shape(exir) = ', shape(exir)
  
  ! optimised effective potential exchange potential
  IF(xctype(1) < 0) CALL oepmain()
  
  ! remove the source term of the exchange-correlation magnetic field if required
  IF(spinpol .and. nosource) CALL projsbf()
  
  ! effective potential from sum of Coulomb and exchange-correlation potentials
  DO ias=1,natmtot
    is=idxis(ias)
    np=npmt(is)
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

include 'my_potcoul.f90'
include 'my_zpotclmt.f90'
include 'my_genzvclmt.f90'
include 'my_zpotcoul.f90'
include 'my_potxc.f90'
include 'my_potxcmt.f90'
include 'my_potxcir.f90'

include 'r_to_zf_mt.f90'
include 'r_to_zf_lm.f90'
include 'z_to_rf_mt.f90'
include 'z_to_rf_lm.f90'
!include 'rf_mt_c_to_f.f90'
!include 'rf_interp.f90'