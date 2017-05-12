SUBROUTINE init_V_ps_loc_G()
  
  USE m_PsPot, ONLY : Ps => Ps_HGH_Params

  USE m_Ps_HGH, ONLY : hgh_eval_Vloc_G

  USE m_hamiltonian, ONLY : V_ps_loc

  USE m_atoms, ONLY : strf => StructureFactor, Nspecies

  USE m_realspace, ONLY : Npoints, Ns
  USE m_PWGrid, ONLY : Gv2
  USE m_cell, ONLY : CellVolume

  IMPLICIT NONE 
  INTEGER :: ip, isp, Nx, Ny, Nz
  REAL(8) :: Gm
  COMPLEX(8), ALLOCATABLE :: ctmp(:)

  WRITE(*,*)
  WRITE(*,*) 'V_ps_loc is initialized via G-space'

  Nx = Ns(1)
  Ny = Ns(2)
  Nz = Ns(3)

  ALLOCATE( ctmp(Npoints) )

  V_ps_loc(:) = 0.d0

  DO isp = 1,Nspecies

    ctmp(:) = cmplx(0.d0,0.d0,kind=8)
    DO ip = 1,Npoints
      Gm = sqrt(Gv2(ip))
      ctmp(ip) = hgh_eval_Vloc_G( Ps(isp), Gm ) * strf(ip,isp) / CellVolume
    ENDDO

    ! inverse FFT: G -> R
    CALL fft_fftw3( ctmp, Nx, Ny, Nz, .true. )

    ! XXX: Move this outside isp loop ?
    DO ip = 1,Npoints
      V_ps_loc(ip) = V_ps_loc(ip) + real( ctmp(ip), kind=8 )
    ENDDO 

  ENDDO 

  DEALLOCATE( ctmp )

END SUBROUTINE 

