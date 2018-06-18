PROGRAM test_PWGrid
  IMPLICIT NONE 
  REAL(8) :: LL(3,3)
  REAL(8) :: ecutwfc_Ry

  ecutwfc_Ry = 30.d0
  LL(1,:) = (/ 16.d0, 0.d0, 0.d0 /)
  LL(2,:) = (/ 0.d0, 16.d0, 0.d0 /)
  LL(3,:) = (/ 0.d0, 0.d0, 16.d0 /)

  CALL init_PWGrid( 0.5d0*ecutwfc_Ry, LL )
  CALL info_PWGrid()

  CALL init_rgrid()

  CALL dealloc_realspace()
  CALL dealloc_PWGrid()

END PROGRAM 

