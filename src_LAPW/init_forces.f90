SUBROUTINE init_forces()
  USE modmain
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

