SUBROUTINE default_kpoints()
  USE m_kpoints, ONLY: reducek, radkpt, ndspem, deltaem, autokpt, vklem, vkloff, ngridk
  IMPLICIT NONE 

  ngridk(:)=1
  vkloff(:)=0.d0
  autokpt=.false.
  radkpt=40.d0
  reducek=1

  vklem(:)=0.d0
  deltaem=0.025d0
  ndspem=1

END SUBROUTINE 