SUBROUTINE default_hamiltonian()
  USE m_hamiltonian, ONLY: tefvr, tefvit, minitefv, maxitefv, epsefvit, befvit, evtype
  IMPLICIT NONE 

  tefvr=.true.
  tefvit=.false.
  minitefv=6
  maxitefv=4
  befvit=0.25d0
  epsefvit=1.d-5
  evtype=1

END SUBROUTINE 