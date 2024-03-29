SUBROUTINE init_forces()
  USE m_atoms, ONLY: natmtot
  USE m_force_stress, ONLY: forcehf, forceibs, forcetot
  IMPLICIT NONE 
  !-------------------------!
  !     force variables     !
  !-------------------------!
  IF( allocated(forcehf) ) DEALLOCATE(forcehf)
  ALLOCATE( forcehf(3,natmtot) )
  IF( allocated(forceibs) ) DEALLOCATE(forceibs)
  ALLOCATE( forceibs(3,natmtot) )
  IF( allocated(forcetot) ) DEALLOCATE(forcetot)
  ALLOCATE( forcetot(3,natmtot) )
END SUBROUTINE

