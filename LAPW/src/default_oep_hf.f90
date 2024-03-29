SUBROUTINE default_oep_hf()
  USE m_oep_hf, ONLY: maxitoep, hybridc, hybrid
  maxitoep=200
  hybrid=.false.
  hybridc=1.d0
END SUBROUTINE