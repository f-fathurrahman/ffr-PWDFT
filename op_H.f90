SUBROUTINE op_H( Ncols, psi, Hpsi )
  USE m_PWGrid, ONLY : Ngwx
  USE m_hamiltonian, ONLY : V_Hartree, V_ps_loc, V_xc
  IMPLICIT NONE 
  INTEGER :: Ncols
  COMPLEX(8) :: psi(Ngwx,Ncols)
  COMPLEX(8) :: Hpsi(Ngwx,Ncols)
  !
  COMPLEX(8), ALLOCATABLE :: Kpsi(:,:), Vpsi(:,:)
  INTEGER :: igw, ic

  ALLOCATE( Kpsi(Ngwx,Ncols) )
  ALLOCATE( Vpsi(Ngwx,Ncols) )

  CALL op_K( Ncols, psi, Kpsi )
  CALL op_V_loc( Ncols, V_Hartree(:) + V_ps_loc(:) + V_xc(:), psi, Vpsi )

  DO ic = 1,Ncols
    DO igw = 1,Ngwx
      Hpsi(igw,ic) = Kpsi(igw,ic) + Vpsi(igw,ic)
    ENDDO 
  ENDDO 

  DEALLOCATE( Kpsi )
  DEALLOCATE( Vpsi )

END SUBROUTINE 


SUBROUTINE op_H_1col( ic, psi, Hpsi )
  USE m_PWGrid, ONLY : Ngwx
  USE m_hamiltonian, ONLY : V_Hartree, V_ps_loc, V_xc
  IMPLICIT NONE 
  COMPLEX(8) :: psi(Ngwx)
  COMPLEX(8) :: Hpsi(Ngwx)
  INTEGER :: ic
  !
  COMPLEX(8), ALLOCATABLE :: Kpsi(:), Vpsi(:)
  INTEGER :: igw
  
  ALLOCATE( Kpsi(Ngwx) )
  ALLOCATE( Vpsi(Ngwx) )

  CALL op_K_1col( psi, Kpsi )
  CALL op_V_loc_1col( V_Hartree(:) + V_ps_loc(:) + V_xc(:), psi, Vpsi )

  DO igw = 1,Ngwx
    Hpsi(igw) = Kpsi(igw) + Vpsi(igw)
  ENDDO 

  DEALLOCATE( Kpsi )
  DEALLOCATE( Vpsi )

END SUBROUTINE 
