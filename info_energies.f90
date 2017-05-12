SUBROUTINE info_energies()
  
  USE m_energies
  IMPLICIT NONE 

  WRITE(*,*)
  WRITE(*,*) 'Energy components:'
  WRITE(*,*)
  WRITE(*,fmt=999) 'Kinetic  = ', E_kinetic
  WRITE(*,fmt=999) 'Local PS = ', E_ps_loc
  WRITE(*,fmt=999) 'NL PS    = ', E_ps_NL
  WRITE(*,fmt=999) 'Hartree  = ', E_Hartree
  WRITE(*,fmt=999) 'XC       = ', E_xc
  WRITE(*,fmt=999) 'NN       = ', E_nn
  WRITE(*,*)       '-----------------------------'
  WRITE(*,fmt=999) 'Total    = ', E_total


  999 FORMAT(1x,A,F18.10)

END SUBROUTINE 
