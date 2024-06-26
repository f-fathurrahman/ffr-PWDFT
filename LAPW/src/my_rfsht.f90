! forward spherical harmonic transform
! convert ... to spherical harmonics
SUBROUTINE my_rfsht(nr,nri,rfmt1,rfmt2)
  USE m_muffin_tins, ONLY: lmmaxi, lmmaxo
  USE m_sht, ONLY: rfshti, rfshto
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr,nri
  REAL(8), intent(in) :: rfmt1(*)
  REAL(8), intent(out) :: rfmt2(*)
  ! local variables
  INTEGER :: nro,i
  
  ! transform the inner part of the muffin-tin
  CALL dgemm('N','N',lmmaxi,nri,lmmaxi,1.d0,rfshti,lmmaxi,rfmt1,lmmaxi,0.d0, &
   rfmt2,lmmaxi)
  
  ! transform the outer part of the muffin-tin
  nro = nr - nri
  i = lmmaxi*nri + 1
  write(*,*) 'my_rfsht: i = ', i
  CALL dgemm('N','N',lmmaxo,nro,lmmaxo,1.d0,rfshto,lmmaxo,rfmt1(i),lmmaxo,0.d0, &
   rfmt2(i),lmmaxo)
  RETURN 
END SUBROUTINE 

