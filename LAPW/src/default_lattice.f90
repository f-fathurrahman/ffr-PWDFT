SUBROUTINE default_lattice()
  USE m_lattice, ONLY: avec, rndavec, epslat
  IMPLICIT NONE 

  avec(:,:)=0.d0
  avec(1,1)=1.d0
  avec(2,2)=1.d0
  avec(3,3)=1.d0
  rndavec=0.d0
  epslat=1.d-6

END SUBROUTINE