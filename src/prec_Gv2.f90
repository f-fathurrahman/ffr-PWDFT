SUBROUTINE prec_Gv2( Ncols, v, Pv )
  USE m_PWGrid, ONLY : Ngwx, idx_gw2r, Gv2
  IMPLICIT NONE
  INTEGER :: Ncols
  COMPLEX(8) :: v(Ngwx,Ncols)
  COMPLEX(8) :: Pv(Ngwx,Ncols)
  INTEGER :: igw, idxG, ic

  DO ic = 1,Ncols
    DO igw = 1,Ngwx
      idxG = idx_gw2r(igw)
      Pv(igw,ic) = v(igw,ic) / (1.d0 + Gv2(idxG))
    ENDDO
  ENDDO

END SUBROUTINE

SUBROUTINE prec_Gv2_inplace( Ncols, v )
  USE m_PWGrid, ONLY : Ngwx, idx_gw2r, Gv2
  IMPLICIT NONE
  INTEGER :: Ncols
  COMPLEX(8) :: v(Ngwx,Ncols)
  INTEGER :: igw, idxG, ic

  DO ic = 1,Ncols
    DO igw = 1,Ngwx
      idxG = idx_gw2r(igw)
      v(igw,ic) = v(igw,ic) / (1.d0 + Gv2(idxG))
    ENDDO
  ENDDO

END SUBROUTINE
