SUBROUTINE default_dos_optics_response()
  USE m_dos_optics_response, ONLY: &
               wsfac, sqados, wplot, nwplot, nswplot, noptcomp, ngrkf, &
               lmirep, lmaxdos, intraband, emaxelnes, dosssum, dosocc, dosmsum, &
               optcomp, vecql
  IMPLICIT NONE 

  sqados(1:2)=0.d0
  sqados(3)=1.d0
  wsfac(1)=-1.d6; wsfac(2)=1.d6
  nwplot=500
  wplot(1)=-0.5d0
  wplot(2)=0.5d0
  nswplot=1

  dosocc=.false.
  dosmsum=.false.
  dosssum=.false.
  lmirep=.true.
  lmaxdos=3
  intraband = .false.
  emaxelnes = -1.2d0
  vecql(:)=0.d0
  noptcomp=1
  optcomp(:,1)=1
  ngrkf=100


END SUBROUTINE