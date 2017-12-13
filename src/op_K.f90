SUBROUTINE op_K( Ncols, psi, Kpsi )
  USE m_PWGrid, ONLY : Ngwx, Gv2, idx_gw2r
  IMPLICIT NONE 
  INTEGER :: Ncols
  COMPLEX(8) :: psi(Ngwx,Ncols)
  COMPLEX(8) :: Kpsi(Ngwx,Ncols)
  INTEGER :: ic, igw, idxG

  DO ic = 1,Ncols
    DO igw = 1,Ngwx
      idxG = idx_gw2r(igw)
      Kpsi(igw,ic) = 0.5d0*Gv2(idxG)*psi(igw,ic)
    ENDDO 
  ENDDO 

END SUBROUTINE 


SUBROUTINE op_K_1col( psi, Kpsi )
  USE m_PWGrid, ONLY : Ngwx, Gv2, idx_gw2r
  IMPLICIT NONE 
  COMPLEX(8) :: psi(Ngwx)
  COMPLEX(8) :: Kpsi(Ngwx)
  INTEGER :: igw, idxG

  DO igw = 1,Ngwx
    idxG = idx_gw2r(igw)
    Kpsi(igw) = 0.5d0*Gv2(idxG)*psi(igw)
  ENDDO 

END SUBROUTINE 


