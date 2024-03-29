SUBROUTINE default_muffin_tins()
  USE m_muffin_tins, ONLY: rmtdelta, rmtall, nrmtscf, lradstp, lmaxo, lmaxi, lmaxapw, &
                     fracinr
  IMPLICIT NONE 

  rmtdelta=0.05d0
  rmtall=-1.d0
  nrmtscf=1.d0
  lradstp=4
  lmaxapw=8
  lmaxo=6
  lmaxi=1
  fracinr=0.01d0

END SUBROUTINE 