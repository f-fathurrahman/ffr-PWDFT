SUBROUTINE my_rotrfmt(rot, nr, nri, rfmt1, rfmt2)
  ! rfmt1 and rfmt2 are packed muffin tin arrays for one atom
  ! rot is a rotation matrix
  use modmain
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(in) :: rot(3,3)
  INTEGER, intent(in) :: nr,nri
  REAL(8), intent(in) :: rfmt1(*)
  REAL(8), intent(out) :: rfmt2(*)
  ! local variables
  INTEGER nro,i
  write(*,*) 'Hello from rotrfmt !!!'
  !
  ! inner part of muffin-tin
  CALL my_rotrflm(rot, lmaxi, nri, lmmaxi, rfmt1, rfmt2)
  ! leading dimension of rfmt1 and rfmt2 is lmmaxi
  !
  ! outer part of muffin-tin
  nro = nr - nri
  i = lmmaxi*nri + 1
  CALL my_rotrflm(rot, lmaxo, nro, lmmaxo, rfmt1(i), rfmt2(i))
  ! leading dimension of rfmt1 and rfmt2 is lmmaxo
  ! 
  RETURN 
END SUBROUTINE 

