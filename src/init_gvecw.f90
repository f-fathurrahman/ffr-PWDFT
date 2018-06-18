SUBROUTINE init_gvecw()
  USE m_PWGrid, ONLY : Ngwx, ecutwfc, Ng, Gv2, idx_gw2r, idx_g2r, Gv, &
                       idx_gw2g
  IMPLICIT NONE 
  INTEGER :: ig, igw
  REAL(8) :: Gvv2

  ! Calculate Ngwx, this should be changed when using non-gamma kpoints
  Ngwx = 0
  DO ig = 1,Ng
    Gvv2 = Gv(1,ig)**2 + Gv(2,ig)**2 + Gv(3,ig)**2
    IF( 0.5d0*Gvv2 <= ecutwfc ) THEN 
      Ngwx = Ngwx + 1
    ENDIF 
  ENDDO 

  ALLOCATE( idx_gw2r(Ngwx) )
  ALLOCATE( idx_gw2g(Ngwx) )
  igw = 0
  DO ig = 1,Ng
    IF( 0.5d0*Gv2(ig) <= ecutwfc ) THEN 
      igw = igw + 1
      idx_gw2r(igw) = idx_g2r(ig)
      idx_gw2g(igw) = ig
    ENDIF 
  ENDDO 

END SUBROUTINE 

