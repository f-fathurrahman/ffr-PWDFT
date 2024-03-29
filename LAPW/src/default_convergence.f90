SUBROUTINE default_convergence()
  USE m_convergence, ONLY: maxscl, epsstress, epspot, epsforce, epsengy
  IMPLICIT NONE 
  maxscl=200
  epspot=1.d-6
  epsengy=1.d-4
  epsforce=5.d-3
  epsstress=1.d-3
END SUBROUTINE