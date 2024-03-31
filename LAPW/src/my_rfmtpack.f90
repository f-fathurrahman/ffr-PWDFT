SUBROUTINE my_rfmtpack(tpack, nr, nri, rfmt1, rfmt2)
  USE m_muffin_tins, ONLY: lmmaxi, lmmaxo
  IMPLICIT NONE 
  ! arguments
  LOGICAL, intent(in) :: tpack
  INTEGER, intent(in) :: nr,nri
  REAL(8), intent(in) :: rfmt1(*)
  REAL(8), intent(out) :: rfmt2(*)
  ! local variables
  INTEGER :: ir,i,j,k
  i = 1
  j = 1

  write(*,*) 'Enter my_rfmtpack'
  write(*,*) 'tpack = ', tpack
  write(*,*) 'nri   = ', nri

  IF(tpack) THEN 
    DO ir = 1,nri
      CALL dcopy(lmmaxi, rfmt1(i), 1, rfmt2(j), 1)
      i = i + lmmaxo
      j = j + lmmaxi
    ENDDO 
  ELSE 
    ! unpacking
    ! zeros rf are added
    DO ir=1,nri
      CALL dcopy(lmmaxi, rfmt1(i), 1, rfmt2(j), 1)
      i = i + lmmaxi
      k = j + lmmaxi
      j = j + lmmaxo
      rfmt2(k:j-1) = 0.d0
    ENDDO 
  ENDIF 
  k = lmmaxo*(nr-nri)
  CALL dcopy(k, rfmt1(i), 1, rfmt2(j), 1)
  RETURN 
END SUBROUTINE 

