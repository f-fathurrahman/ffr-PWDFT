SUBROUTINE prec_Gv2( Ncols, v, Pv )
  USE m_PWGrid, ONLY : Ngwx, idx_gw2g, Gv2
  IMPLICIT NONE
  INTEGER :: Ncols
  COMPLEX(8) :: v(Ngwx,Ncols)
  COMPLEX(8) :: Pv(Ngwx,Ncols)
  INTEGER :: igw, ig, ic

  DO ic = 1,Ncols
    DO igw = 1,Ngwx
      ig = idx_gw2g(igw)
      Pv(igw,ic) = v(igw,ic) / (1.d0 + Gv2(ig))
    ENDDO
  ENDDO

END SUBROUTINE

SUBROUTINE prec_Gv2_inplace( Ncols, v )
  USE m_PWGrid, ONLY : Ngwx, idx_gw2g, Gv2
  IMPLICIT NONE
  INTEGER :: Ncols
  COMPLEX(8) :: v(Ngwx,Ncols)
  INTEGER :: igw, ig, ic

  DO ic = 1,Ncols
    DO igw = 1,Ngwx
      ig = idx_gw2g(igw)
      v(igw,ic) = v(igw,ic) / (1.d0 + Gv2(ig))
    ENDDO
  ENDDO

END SUBROUTINE
