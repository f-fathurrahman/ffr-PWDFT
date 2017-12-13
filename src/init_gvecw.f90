SUBROUTINE init_gvecw()
  USE m_PWGrid, ONLY : Ngwx, ecutwfc, Ng, Gv2, idx_gw2r
  IMPLICIT NONE 
  INTEGER :: ig, igw

  Ngwx = 0
  DO ig = 1,Ng
    IF( 0.5d0*Gv2(ig) <= ecutwfc ) THEN 
      Ngwx = Ngwx + 1
    ENDIF 
  ENDDO 

  ALLOCATE( idx_gw2r(Ngwx) )
  igw = 0
  DO ig = 1,Ng
    IF( 0.5d0*Gv2(ig) <= ecutwfc ) THEN 
      igw = igw + 1
      idx_gw2r(igw) = ig
    ENDIF 
  ENDDO 


END SUBROUTINE 

