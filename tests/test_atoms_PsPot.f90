PROGRAM test_atoms_PsPot

  USE m_PsPot, ONLY : PsPot_Dir
  IMPLICIT NONE 
  INTEGER :: Narg
  CHARACTER(64) :: filexyz

  Narg = iargc()
  IF( Narg /= 1 ) THEN 
    WRITE(*,*) 'ERROR: exactly one argument must be given'
    STOP 
  ENDIF 

  CALL getarg(1,filexyz)

  CALL init_atoms_xyz(filexyz)

  ! Override PsPot_Dir
  PsPot_Dir = '../HGH/'
  CALL init_PsPot()

  CALL info_atoms()
  CALL info_PsPot()

  CALL dealloc_PsPot()
  CALL dealloc_atoms()

END PROGRAM 
