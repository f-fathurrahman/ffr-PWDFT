SUBROUTINE default_apwlo()
  USE m_apwlo, ONLY: nxoapwlo, lorbcnd, lorbordc, nxlo, epsband, dlefe, demaxbnd, deapwlo, &
               autolinengy
  IMPLICIT NONE 

  nxoapwlo=0
  nxlo=0
  lorbcnd=.false.
  lorbordc=3
  epsband=1.d-12
  demaxbnd=2.5d0
  autolinengy=.false.
  dlefe=-0.1d0
  deapwlo=0.05d0

END SUBROUTINE