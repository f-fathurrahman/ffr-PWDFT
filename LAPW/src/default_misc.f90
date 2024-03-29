SUBROUTINE default_misc()
  USE m_misc, ONLY: wrtvars, scrpath, nwrite, ntasks, notelns
  IMPLICIT NONE 

  wrtvars=.false.
  scrpath=''
  nwrite=0
  ntasks=0
  notelns=0

END SUBROUTINE