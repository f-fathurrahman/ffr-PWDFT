!-------------------------------
SUBROUTINE init_charges_states()
!-------------------------------
  USE m_constants, ONLY: fourpi
  USE m_atoms, ONLY: nspecies, spzn, natoms, nstsp, nstspmax, occsp, ksp, spze, spcore, nstcr
  USE m_lattice, ONLY: omega
  USE m_charge_moment_current, ONLY: chgval, chgtot, chgzn, chgexs, chgcrtot, chgcr, rwigner
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

END SUBROUTINE 