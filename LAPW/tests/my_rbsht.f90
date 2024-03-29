SUBROUTINE my_rbsht(nr,nri,rfmt1,rfmt2)
  USE m_muffin_tins, ONLY: lmmaxi, lmmaxo
  USE m_sht, ONLY: rbshti, rbshto
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr,nri
  REAL(8), intent(in) :: rfmt1(*)
  REAL(8), intent(out) :: rfmt2(*)
  ! local variables
  INTEGER :: nro,i

  ! transform the inner part of the muffin-tin

  write(*,*)
  write(*,*) 'my_rbsht: Inner part of MT'
  write(*,*) 'my_rbsht: Matrix A dimension = ', lmmaxi, lmmaxi
  write(*,*) 'my_rbsht: Matrix B dimension = ', lmmaxi, nri
  write(*,*) 'my_rbsht: Matrix C dimension = ', lmmaxi, nri
  write(*,*) 'Number of elements: lmmaxi*nri = ', lmmaxi*nri
  CALL dgemm('N','N',lmmaxi,nri,lmmaxi,1.d0,rbshti,lmmaxi,rfmt1,lmmaxi,0.d0, &
   rfmt2,lmmaxi)

  write(102,*) rfmt1(1:lmmaxi*nri)
  write(103,*) rfmt2(1:lmmaxi*nri)
  
  ! transform the outer part of the muffin-tin
  nro=nr-nri
  i=lmmaxi*nri+1
  write(*,*)
  write(*,*) 'my_rbsht: Outer part of MT'
  write(*,*) 'my_rbsht: i   = ', i
  write(*,*) 'my_rbsht: Matrix A dimension = ', lmmaxo, lmmaxo
  write(*,*) 'my_rbsht: Matrix B dimension = ', lmmaxo, nro
  write(*,*) 'my_rbsht: Matrix C dimension = ', lmmaxo, nro
  write(*,*) 'my_rbsht: Number of elements: lmmaxo*nro = ', lmmaxo*nro
  write(*,*) 'my_rbsht: Last index accessed = ', i + lmmaxo*nro
  CALL dgemm('N','N',lmmaxo,nro,lmmaxo,1.d0,rbshto,lmmaxo,rfmt1(i),lmmaxo,0.d0, &
   rfmt2(i),lmmaxo)

  write(104,*) rfmt1(i:i+lmmaxo*nro-1)
  write(105,*) rfmt2(i:i+lmmaxo*nro-1)

  !call write_radial_mt()
  !stop 'should plot rfmt1 and rfmt2'

  RETURN 
END SUBROUTINE 

