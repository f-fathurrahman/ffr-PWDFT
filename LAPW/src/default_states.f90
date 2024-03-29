SUBROUTINE default_states()
  USE m_states, ONLY: swidth, stype, scissor, nkstlist, nempty0, mstar, evaltol, &
                epsocc, autoswidth, kstlist
  IMPLICIT NONE 

  swidth = 0.001d0
  autoswidth = .false.
  stype = 3
  scissor = 0.d0
  nkstlist = 1
  kstlist(:,1) = 1
  nempty0 = 4.d0
  mstar = 10.d0
  evaltol = -1.d0
  epsocc = 1.d-8

END SUBROUTINE 