!-------------------------------
SUBROUTINE init_charges_states()
!-------------------------------
  
  USE modmain, ONLY: &
               nspecies, spzn, natoms, nstsp, nstspmax, chgcr, occsp, &
               ksp, spcore, spze, rwigner, omega, nstcr, fourpi, &
               chgval, chgtot, chgzn, chgexs, chgcrtot

  IMPLICIT NONE 
  INTEGER :: is, ist

  !--------------------------------------!
  !     charges and number of states     !
  !--------------------------------------!
  chgzn = 0.d0
  chgcrtot = 0.d0
  chgval = 0.d0
  nstspmax = 0
  nstcr = 0
  DO is = 1,nspecies
    ! nuclear charge
    chgzn = chgzn + spzn(is)*natoms(is)
    ! find the maximum number of atomic states
    nstspmax = max(nstspmax,nstsp(is))
    ! compute the electronic charge for each species, as well as the total core and
    ! valence charge
    spze(is) = 0.d0
    chgcr(is) = 0.d0
    DO ist = 1,nstsp(is)
      spze(is) = spze(is) + occsp(ist,is)
      IF( spcore(ist,is) ) THEN 
        chgcr(is) = chgcr(is) + occsp(ist,is)
        nstcr = nstcr + 2*ksp(ist,is)*natoms(is)
      ELSE 
        chgval = chgval + occsp(ist,is)*natoms(is)
      ENDIF 
    ENDDO 
    chgcrtot=chgcrtot + chgcr(is)*natoms(is)
  ENDDO 
  ! add excess charge
  chgval = chgval + chgexs
  ! total charge
  chgtot = chgcrtot + chgval
  IF( chgtot < 1.d-8 ) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(init0): zero total charge")')
    WRITE(*,*)
    STOP 
  ENDIF 
  
  ! effective Wigner radius
  rwigner = (3.d0/(fourpi*(chgtot/omega)))**(1.d0/3.d0)
  
  ! write to VARIABLES.OUT
  !call writevars('spze',nv=nspecies,rva=spze)
  !call writevars('chgcr',nv=nspecies,rva=chgcr)
  !call writevars('chgexs',rv=chgexs)
  !call writevars('chgval',rv=chgtot)

END SUBROUTINE 

