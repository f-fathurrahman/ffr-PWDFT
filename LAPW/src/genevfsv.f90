SUBROUTINE genevfsv()
  USE m_states, ONLY: nstsv, nstfv, evalsv
  USE m_spin, ONLY: nspnfv
  USE m_hamiltonian, ONLY: nmatmax
  USE m_kpoints, ONLY: nkpt
  USE m_misc, ONLY: filext
  IMPLICIT NONE 
  ! local variables
  INTEGER :: ik,lp
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: evalfv(:,:)
  COMPLEX(8), ALLOCATABLE :: evecfv(:,:,:),evecsv(:,:)
  
  ALLOCATE(evalfv(nstfv,nspnfv))
  ALLOCATE(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
  
  DO ik=1,nkpt
    ! solve the first- and second-variational eigenvalue equations
    CALL eveqn(ik,evalfv,evecfv,evecsv)
    ! write the eigenvalues/vectors to file
    CALL putevalfv(filext,ik,evalfv)
    CALL putevalsv(filext,ik,evalsv(:,ik))
    CALL putevecfv(filext,ik,evecfv)
    CALL putevecsv(filext,ik,evecsv)
  ENDDO 
  DEALLOCATE(evalfv,evecfv,evecsv)

  RETURN 
END SUBROUTINE 

