SUBROUTINE op_V_loc( Ncols, V_loc, psi, Vpsi )
  USE m_realspace, ONLY : Npoints, Ns
  USE m_PWGrid, ONLY : Ngwx, idx_gw2r
  IMPLICIT NONE 
  INTEGER :: Ncols
  REAL(8) :: V_loc(Npoints)
  COMPLEX(8) :: psi(Ngwx,Ncols)
  COMPLEX(8) :: Vpsi(Ngwx,Ncols)
  !
  COMPLEX(8), ALLOCATABLE :: ctmp(:)
  INTEGER :: igw, ip, ic

  ALLOCATE( ctmp(Npoints) )

  DO ic = 1,Ncols
    ctmp(:) = cmplx(0.d0,0.d0,kind=8)
    DO igw = 1,Ngwx
      ip = idx_gw2r(igw)
      ctmp(ip) = psi(igw,ic)
    ENDDO 
    ! transform to real space via inverse FFT
    CALL fft_fftw3( ctmp, Ns(1), Ns(2), Ns(3), .TRUE. )
    !
    ! Multiplication in real space
    DO ip = 1,Npoints
      ctmp(ip) = ctmp(ip)*V_loc(ip)
    ENDDO 
    !
    CALL fft_fftw3( ctmp, Ns(1), Ns(2), Ns(3), .FALSE. )
    !
    DO igw = 1,Ngwx
      ip = idx_gw2r(igw)
      Vpsi(igw,ic) = ctmp(ip)
    ENDDO 
  ENDDO 

END SUBROUTINE 


SUBROUTINE op_V_loc_1col( V_loc, psi, Vpsi )
  USE m_realspace, ONLY : Npoints, Ns
  USE m_PWGrid, ONLY : Ngwx, idx_gw2r
  IMPLICIT NONE 
  REAL(8) :: V_loc(Npoints)
  COMPLEX(8) :: psi(Ngwx)
  COMPLEX(8) :: Vpsi(Ngwx)
  !
  COMPLEX(8), ALLOCATABLE :: ctmp(:)
  INTEGER :: igw, ip

  ALLOCATE( ctmp(Npoints) )

  ctmp(:) = cmplx(0.d0,0.d0,kind=8)
  DO igw = 1,Ngwx
    ip = idx_gw2r(igw)
    ctmp(ip) = psi(igw)
  ENDDO 
  ! transform to real space via inverse FFT
  CALL fft_fftw3( ctmp, Ns(1), Ns(2), Ns(3), .TRUE. )
  !
  ! Multiplication in real space
  DO ip = 1,Npoints
    ctmp(ip) = ctmp(ip)*V_loc(ip)
  ENDDO 
  !
  CALL fft_fftw3( ctmp, Ns(1), Ns(2), Ns(3), .FALSE. )
  !
  DO igw = 1,Ngwx
    ip = idx_gw2r(igw)
    Vpsi(igw) = ctmp(ip)
  ENDDO 

END SUBROUTINE 

