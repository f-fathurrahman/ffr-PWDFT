SUBROUTINE rfmtlm(lm,nr,nri,rfmt,fr)
  USE m_muffin_tins, ONLY: nrmtmax, npmtmax, lmmaxo, lmmaxi
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: lm,nr,nri
  REAL(8), intent(in) :: rfmt(npmtmax)
  REAL(8), intent(out) :: fr(nrmtmax)
  ! local variables
  INTEGER iro,ir,i

  IF(lm > lmmaxi) THEN 
    fr(1:nri)=0.d0
  ELSE 
    i=lm
    DO ir=1,nri
      fr(ir)=rfmt(i)
      i=i+lmmaxi
    ENDDO 
  ENDIF 
  
  iro=nri+1
  IF(lm > lmmaxo) THEN 
    fr(iro:nr)=0.d0
  ELSE 
    i=lmmaxi*nri+lm
    DO ir=iro,nr
      fr(ir)=rfmt(i)
      i=i+lmmaxo
    ENDDO 
  ENDIF 

  RETURN 
END SUBROUTINE 

