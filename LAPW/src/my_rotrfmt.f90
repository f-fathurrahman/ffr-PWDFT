SUBROUTINE my_rotrfmt(rot, nr, nri, rfmt1, rfmt2)
  use modmain
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(in) :: rot(3,3)
  INTEGER, intent(in) :: nr,nri
  REAL(8), intent(in) :: rfmt1(*)
  REAL(8), intent(out) :: rfmt2(*)
  ! local variables
  INTEGER nro,i
  ! inner part of muffin-tin
  CALL my_rotrflm(rot, lmaxi, nri, lmmaxi, rfmt1, rfmt2)
  ! outer part of muffin-tin
  nro = nr-nri
  i = lmmaxi*nri + 1
  CALL my_rotrflm(rot, lmaxo, nro, lmmaxo, rfmt1(i), rfmt2(i))
  ! what is the leading dimension of rfmt1 and rftm2?
  RETURN 
END SUBROUTINE 

