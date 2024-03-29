!----------------------------
SUBROUTINE init_atoms_cores()
!----------------------------
  USE m_constants, ONLY: y00
  USE m_atoms, ONLY: nspecies, natmtot, idxis
  USE m_atomic_species, ONLY: ptnucl, rsp, spzn, vcln, nstspmax, nstsp, evalsp, occsp, &
               nrsp, nrspmax
  USE m_core_states, ONLY: spincore, evalcr, rhocr, occcr, nspncr, rwfcr
  USE m_muffin_tins, ONLY: nrmtmax
  IMPLICIT NONE 
  INTEGER :: is, nr, ist, ias
  REAL(8) :: t1

  !-------------------------!
  !     atoms and cores     !
  !-------------------------!
  
  ! determine the nuclear Coulomb potential
  IF( allocated(vcln) ) DEALLOCATE(vcln)
  ALLOCATE( vcln(nrspmax,nspecies) )
  t1 = 1.d0/y00
  DO is = 1,nspecies
    nr = nrsp(is)
    CALL potnucl( ptnucl, nr, rsp(:,is), spzn(is), vcln(:,is) )
    vcln(1:nr,is) = t1*vcln(1:nr,is)
  ENDDO 
  
  ! solve the Kohn-Sham-Dirac equations for all atoms
  CALL allatoms()
  
  ! allocate core state occupancy and eigenvalue arrays and set to default
  IF( allocated(occcr)) DEALLOCATE(occcr)
  ALLOCATE( occcr(nstspmax, natmtot) )
  IF( allocated(evalcr) ) DEALLOCATE(evalcr)
  ALLOCATE( evalcr(nstspmax, natmtot) )
  DO ias = 1,natmtot
    is = idxis(ias)
    DO ist = 1,nstsp(is)
      occcr(ist,ias) = occsp(ist,is)
      evalcr(ist,ias) = evalsp(ist,is)
    ENDDO 
  ENDDO 
  
  ! allocate core state radial wavefunction array
  IF( allocated(rwfcr) ) DEALLOCATE(rwfcr)
  ALLOCATE( rwfcr(nrspmax, 2, nstspmax, natmtot) )
  
  ! number of core spin channels
  IF( spincore ) THEN 
    nspncr = 2
  ELSE 
    nspncr = 1
  ENDIF 
  
  ! allocate core state charge density array
  IF( allocated(rhocr) ) DEALLOCATE(rhocr)
  ALLOCATE( rhocr(nrmtmax, natmtot, nspncr) )

END SUBROUTINE 

