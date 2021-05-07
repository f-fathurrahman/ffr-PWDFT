SUBROUTINE default_density_pot_xc()
  USE m_density_pot_xc, ONLY: xctype, trhonorm, tc_tb09, taudft, nosource, msmooth, maxitksi, &
                dncgga, c_tb09, tauksi
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

END SUBROUTINE 