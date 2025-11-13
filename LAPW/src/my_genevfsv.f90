SUBROUTINE my_genevfsv()
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
  COMPLEX(8), ALLOCATABLE :: evecfv(:,:,:),evecsv(:,:)

  write(*,*)
  write(*,*) '<div> ENTER: my_genevfsv'
  write(*,*)

  ALLOCATE(evalfv(nstfv,nspnfv))
  ALLOCATE(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
  
  DO ik=1,nkpt
    ! solve the first- and second-variational eigenvalue equations
    CALL my_eveqn(ik,evalfv,evecfv,evecsv)
    ! write the eigenvalues/vectors to file
    CALL putevalfv(filext, ik, evalfv)
    CALL putevalsv(filext, ik, evalsv(:,ik))
    CALL putevecfv(filext, ik, evecfv)
    CALL putevecsv(filext, ik, evecsv)
  ENDDO
  ! 
  ! The eigenvalues and eigenvectors are already written to files
  ! Temporary arrays are deallocated here.
  DEALLOCATE(evalfv,evecfv,evecsv)

  write(*,*)
  write(*,*) '</div> EXIT: my_genevfsv'
  write(*,*)

  RETURN 
END SUBROUTINE 



!!-----------------------
!SUBROUTINE my_genevfsv()
!!-----------------------
!  USE m_states, ONLY: nstsv, nstfv, evalsv
!  USE m_spin, ONLY: nspnfv
!  USE m_hamiltonian, ONLY: nmatmax
!  USE m_kpoints, ONLY: nkpt
!  USE m_misc, ONLY: filext
!  IMPLICIT NONE 
!  ! local variables
!  INTEGER :: ik, i
!  ! ALLOCATABLE arrays
!  REAL(8), ALLOCATABLE :: evalfv(:,:)
!  COMPLEX(8), ALLOCATABLE :: evecfv(:,:,:), evecsv(:,:)
!
!  write(*,*)
!  write(*,*) '<div> ENTER: my_genevfsv'
!  write(*,*)
!  write(*,*) 'NOTE:'
!  write(*,*) 'NOTE: Second variational eigenvalues are not calculated'
!  write(*,*) 'NOTE:'
!  write(*,*)
!
!  ! XXXX Only nspnfv == 1 is tested
!
!  ALLOCATE(evalfv(nstfv,nspnfv))
!  ALLOCATE(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
!  evalfv(:,:) = 0.0d0
!  evecfv(:,:,:) = 0.0d0
!  evecsv(:,:) = 0.d0
!
!  DO ik=1,nkpt
!    ! solve the first- and second-variational eigenvalue equations
!    CALL my_eveqn(ik,evalfv,evecfv)
!    ! no calculation of second-variational eigenvectors
!    !
!    !XXX If we don't compute 2nd variational eigenvectors, the evalsv will be
!    !XXX written using evalfv in the following subroutine call.
!    !XXX evecsv is set to diagonal 1.
!    !XXX evalsv will be used in occupy. evalfv is not stored as global variable
!    !
!    DO i = 1,nstsv
!      evalsv(i,ik) = evalfv(i,1)
!    ENDDO 
!    evecsv(:,:) = 0.d0
!    DO i = 1,nstsv
!      evecsv(i,i) = 1.d0
!    ENDDO
!    !
!    ! write the eigenvalues/vectors to file
!    CALL putevalfv(filext, ik, evalfv)
!    CALL putevalsv(filext, ik, evalsv(:,ik))
!    CALL putevecfv(filext, ik, evecfv)
!    CALL putevecsv(filext, ik, evecsv)
!  ENDDO
!  ! 
!  ! The eigenvalues and eigenvectors are already written to files
!  ! Temporary arrays are deallocated here.
!  DEALLOCATE(evalfv,evecfv,evecsv)
!
!  write(*,*)
!  write(*,*) '</div> EXIT: my_genevfsv'
!  write(*,*)
!
!  RETURN 
!END SUBROUTINE 

