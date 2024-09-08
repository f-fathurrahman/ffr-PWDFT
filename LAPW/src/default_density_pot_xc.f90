SUBROUTINE default_density_pot_xc()
  USE m_density_pot_xc
  IMPLICIT NONE 

  xctype(1)=3
  xctype(2:3)=0
  trhonorm=.true.
  c_tb09=0.d0
  tc_tb09=.false.
  maxitksi=200
  taudft=.false.
  tauksi=0.002d0
  msmooth=0
  dncgga=1.d-8
  nosource=.false.
  ssxc=1.d0

END SUBROUTINE 