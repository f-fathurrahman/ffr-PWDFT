SUBROUTINE op_K( Ncols, psi, Kpsi )
  USE m_PWGrid, ONLY : Ngwx, Gv2, idx_gw2g
  IMPLICIT NONE 
  INTEGER :: Ncols
  COMPLEX(8) :: psi(Ngwx,Ncols)
  COMPLEX(8) :: Kpsi(Ngwx,Ncols)
  INTEGER :: ic, igw, ig

  DO ic = 1,Ncols
    DO igw = 1,Ngwx
      ig = idx_gw2g(igw)
      Kpsi(igw,ic) = 0.5d0*Gv2(ig)*psi(igw,ic)
    ENDDO 
  ENDDO 

END SUBROUTINE 


SUBROUTINE op_K_1col( psi, Kpsi )
  USE m_PWGrid, ONLY : Ngwx, Gv2, idx_gw2g
  IMPLICIT NONE 
  COMPLEX(8) :: psi(Ngwx)
  COMPLEX(8) :: Kpsi(Ngwx)
  INTEGER :: igw, ig

  DO igw = 1,Ngwx
    ig = idx_gw2g(igw)
    Kpsi(igw) = 0.5d0*Gv2(ig)*psi(igw)
  ENDDO 

END SUBROUTINE 


