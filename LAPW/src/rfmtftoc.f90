SUBROUTINE rfmtftoc(nrc,nrci,rfmt,rfcmt)
  USE modmain
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nrc,nrci
  REAL(8), intent(in) :: rfmt(*)
  REAL(8), intent(out) :: rfcmt(*)
  ! local variables
  INTEGER irc,i,j,n
  i=1
  j=1
  n=lmmaxi*lradstp
  DO irc=1,nrci
    CALL dcopy(lmmaxi,rfmt(i),1,rfcmt(j),1)
    i=i+n
    j=j+lmmaxi
  ENDDO 
  i=i+(lradstp-1)*(lmmaxo-lmmaxi)
  n=lmmaxo*lradstp
  DO irc=nrci+1,nrc
    CALL dcopy(lmmaxo,rfmt(i),1,rfcmt(j),1)
    i=i+n
    j=j+lmmaxo
  ENDDO 
  RETURN 
END SUBROUTINE 