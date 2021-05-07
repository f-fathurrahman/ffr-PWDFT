!------------------------------
SUBROUTINE init_chgden_pot_xc()
!------------------------------

  USE modmain

  IMPLICIT NONE 

  !-------------------------------------------------------------!
  !     charge density, potentials and exchange-correlation     !
  !-------------------------------------------------------------!
  
  ! allocate charge density arrays
  IF( allocated(rhomt) ) DEALLOCATE(rhomt)
  ALLOCATE( rhomt(npmtmax, natmtot) )
  IF( allocated(rhoir) ) DEALLOCATE(rhoir)
  ALLOCATE( rhoir(ngtot) )
  
  ! allocate magnetisation arrays
  IF( allocated(magmt) ) DEALLOCATE(magmt)
  IF( allocated(magir) ) DEALLOCATE(magir)
  IF( spinpol ) THEN 
    ALLOCATE( magmt(npmtmax, natmtot, ndmag) )
    ALLOCATE( magir(ngtot, ndmag) )
  ENDIF 
  
  ! check if the current density should be calculated
  IF( tafield .or. ( any(task == [371,372,373,460,461])) ) THEN 
    tcden = .true.
  ENDIF 
  
  ! allocate current density arrays
  IF( allocated(cdmt) ) DEALLOCATE(cdmt)
  IF( allocated(cdir) ) DEALLOCATE(cdir)
  IF( tcden ) THEN 
    ALLOCATE(cdmt(npmtmax,natmtot,3),cdir(ngtot,3))
  ENDIF 
  
  ! Coulomb potential
  IF( allocated(vclmt)) DEALLOCATE(vclmt)
  ALLOCATE( vclmt(npmtmax,natmtot) )
  IF( allocated(vclir) ) DEALLOCATE( vclir )
  ALLOCATE( vclir(ngtot) )
  
  ! exchange energy density
  IF( allocated(exmt) ) DEALLOCATE(exmt)
  ALLOCATE( exmt(npmtmax,natmtot) )
  IF( allocated(exir) ) DEALLOCATE(exir)
  ALLOCATE( exir(ngtot) )
  
  ! correlation energy density
  IF( allocated(ecmt) ) DEALLOCATE(ecmt)
  ALLOCATE( ecmt(npmtmax,natmtot) )
  IF( allocated(ecir) ) DEALLOCATE(ecir)
  ALLOCATE( ecir(ngtot) )
  
  ! exchange-correlation potential
  IF( allocated(vxcmt) ) DEALLOCATE(vxcmt)
  ALLOCATE( vxcmt(npmtmax,natmtot) )
  IF( allocated(vxcir) ) DEALLOCATE(vxcir)
  ALLOCATE( vxcir(ngtot) )
  
  ! effective Kohn-Sham potential
  IF( allocated(vsmt) ) DEALLOCATE(vsmt)
  ALLOCATE( vsmt(npmtmax,natmtot) )
  IF( allocated(vsir) ) DEALLOCATE(vsir)
  ALLOCATE(vsir(ngtot))
  IF( allocated(vsig) ) DEALLOCATE(vsig)
  ALLOCATE( vsig(ngvec) )
  
  ! exchange-correlation, dipole and Kohn-Sham effective magnetic fields
  IF( allocated(bxcmt) ) DEALLOCATE(bxcmt)
  IF( allocated(bxcir) ) DEALLOCATE(bxcir)
  IF( allocated(bdmt) ) DEALLOCATE(bdmt)
  IF( allocated(bdir) ) DEALLOCATE(bdir)
  IF( allocated(bsmt) ) DEALLOCATE(bsmt)
  IF( allocated(bsir) ) DEALLOCATE(bsir)
  IF( spinpol ) THEN 
    ALLOCATE( bxcmt(npmtmax,natmtot,ndmag), bxcir(ngtot,ndmag) )
    IF(tbdip) ALLOCATE(bdmt(npmtmax,natmtot,ndmag), bdir(ngtot,ndmag) )
    ALLOCATE( bsmt(npcmtmax,natmtot,ndmag), bsir(ngtot,ndmag) )
  ENDIF 
  
  ! kinetic energy density
  IF( allocated(taumt) ) DEALLOCATE( taumt )
  IF( allocated(tauir) ) DEALLOCATE( tauir )
  IF( allocated(taucr) ) DEALLOCATE( taucr )
  IF( (xcgrad == 3) .or. (xcgrad == 4)) THEN 
    ALLOCATE( taumt(npmtmax,natmtot,nspinor), tauir(ngtot,nspinor) )
    ALLOCATE( taucr(npmtmax,natmtot,nspinor) )
  ENDIF 
  
  ! tau-DFT exchange-correlation and Kohn-Sham potentials
  IF( allocated(wxcmt) ) DEALLOCATE(wxcmt)
  IF( allocated(wxcir) ) DEALLOCATE(wxcir)
  IF( allocated(wsmt) ) DEALLOCATE(wsmt)
  IF( allocated(wsir) ) DEALLOCATE(wsir)
  IF( xcgrad == 4 ) THEN 
    ALLOCATE( wxcmt(npmtmax,natmtot), wxcir(ngtot) )
    ALLOCATE( wsmt(npcmtmax,natmtot), wsir(ngtot) )
    tevecsv = .true.
  ENDIF 
  
  ! spin-orbit coupling radial function
  IF( allocated(socfr) ) DEALLOCATE(socfr)
  IF( spinorb ) THEN 
    ALLOCATE( socfr(nrcmtmax, natmtot) )
  ENDIF 
  
  ! allocate muffin-tin charge and moment arrays
  IF(allocated(chgcrlk) ) DEALLOCATE(chgcrlk)
  ALLOCATE( chgcrlk(natmtot) )
  IF( allocated(chgmt) ) DEALLOCATE(chgmt)
  ALLOCATE( chgmt(natmtot) )
  IF( allocated(mommt) ) DEALLOCATE(mommt)
  ALLOCATE( mommt(3,natmtot) )
  
  ! check if scaled spin exchange-correlation (SSXC) should be used
  IF( abs(ssxc-1.d0) > 1.d-6 ) THEN 
    tssxc = .true.
  ELSE 
    tssxc = .false.
  ENDIF 

END SUBROUTINE 

