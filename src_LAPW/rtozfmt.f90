SUBROUTINE rtozfmt(nr,nri,rfmt,zfmt)
  USE m_muffin_tins, ONLY: lmaxi, lmmaxo, lmmaxi, lmaxo
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr,nri
  REAL(8), intent(in) :: rfmt(*)
  COMPLEX(8), intent(out) :: zfmt(*)
  ! local variables
  INTEGER ir,i
  i=1
  DO ir=1,nri
    CALL rtozflm(lmaxi,rfmt(i),zfmt(i))
    i=i+lmmaxi
  ENDDO 
  DO ir=nri+1,nr
    CALL rtozflm(lmaxo,rfmt(i),zfmt(i))
    i=i+lmmaxo
  ENDDO 
  RETURN 
END SUBROUTINE 

