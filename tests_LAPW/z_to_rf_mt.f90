SUBROUTINE z_to_rf_mt(nr, nri, zfmt, rfmt)
  USE m_muffin_tins, ONLY: lmmaxi, lmmaxo, lmaxo, lmaxi
  IMPLICIT NONE 
  ! arguments
  INTEGER, INTENT(in) :: nr,nri
  COMPLEX(8), INTENT(in) :: zfmt(*)
  REAL(8), INTENT(out) :: rfmt(*)
  ! local variables
  INTEGER :: ir,i

  i = 1
  DO ir = 1,nri
    CALL z_to_rf_lm(lmaxi, zfmt(i), rfmt(i))
    i = i + lmmaxi
  ENDDO 
  DO ir = nri+1,nr
    CALL z_to_rf_lm(lmaxo, zfmt(i), rfmt(i))
    i = i + lmmaxo
  ENDDO 
  RETURN 
END SUBROUTINE 