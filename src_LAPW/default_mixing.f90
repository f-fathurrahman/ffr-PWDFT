SUBROUTINE default_mixing()
  USE m_mixing, ONLY: broydpm, amixpm, mixtype, mixsdb
  
  ! Broyden parameters recommended by M. Meinert
  mixsdb=5
  broydpm(1)=0.4d0
  broydpm(2)=0.15d0

  mixtype=3
  amixpm(1)=0.05d0
  amixpm(2)=1.d0

END SUBROUTINE 