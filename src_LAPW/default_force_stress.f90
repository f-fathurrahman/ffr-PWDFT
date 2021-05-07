SUBROUTINE default_force_stress()
  USE m_force_stress, ONLY: tforce, tfav0, tau0latv, tau0atp, maxlatvstp, maxatpstp, latvopt, &
              deltast, atpopt
  IMPLICIT NONE 
  
  tforce=.false.
  tfav0=.true.
  tau0latv=0.2d0
  tau0atp=0.25d0
  maxlatvstp=30
  deltast=0.001d0
  atpopt=1
  maxatpstp=200
  latvopt=0
END SUBROUTINE 
