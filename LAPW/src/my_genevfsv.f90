!-----------------------
SUBROUTINE my_genevfsv()
!-----------------------
  USE m_states, ONLY: nstsv, nstfv, evalsv
  USE m_spin, ONLY: nspnfv
  USE m_hamiltonian, ONLY: nmatmax
  USE m_kpoints, ONLY: nkpt
  USE m_misc, ONLY: filext
  IMPLICIT NONE 
  ! local variables
  INTEGER :: ik
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: evalfv(:,:)
  COMPLEX(8), ALLOCATABLE :: evecfv(:,:,:), evecsv(:,:)

  write(*,*)
  write(*,*) '****** ENTER: my_genevfsv'
  write(*,*)
  write(*,*) 'NOTE:'
  write(*,*) 'NOTE: Second variational eigenvalues are not calculated'
  write(*,*) 'NOTE:'
  write(*,*)

  ALLOCATE(evalfv(nstfv,nspnfv))
  ALLOCATE(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
  evalfv(:,:) = 0.0d0
  evecfv(:,:,:) = 0.0d0
  evecsv(:,:) = 0.d0

  DO ik=1,nkpt
    ! solve the first- and second-variational eigenvalue equations
    CALL my_eveqn(ik,evalfv,evecfv)
    !
    ! write the eigenvalues/vectors to file
    CALL putevalfv(filext,ik,evalfv)
    CALL putevalsv(filext,ik,evalsv(:,ik))
    CALL putevecfv(filext,ik,evecfv)
    CALL putevecsv(filext,ik,evecsv)
  ENDDO
  ! 
  ! The eigenvalues and eigenvectors are already written to files
  ! Temporary arrays are deallocated here.
  DEALLOCATE(evalfv,evecfv,evecsv)

  write(*,*)
  write(*,*) '****** EXIT: my_genevfsv'
  write(*,*)

  RETURN 
END SUBROUTINE 

