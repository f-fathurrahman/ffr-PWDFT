! Fadjar Fathurrahman (20910015), May 2011
! Updated: February 2017
!
! Calling FFTW from legacy Fortran.
! In place FFT using FFTW3: `zdata` is overwritten
!-----------------------------------------------------
SUBROUTINE fft_fftw3( zdata, Nx, Ny, Nz, t_inv)
!-----------------------------------------------------
  IMPLICIT NONE
  INCLUDE 'fftw3.f'
  ! Arguments
  INTEGER :: Nx, Ny, Nz
  COMPLEX(8) :: zdata(Nx, Ny, Nz)
  LOGICAL :: t_inv

  ! TODO: Probably it is better to save PLAN_BACKWARD and PLAN_FORWARD
  !       as global variables
  INTEGER(8) :: PLAN_BACKWARD
  INTEGER(8) :: PLAN_FORWARD
 
  IF( t_inv ) THEN ! backward transform
    CALL dfftw_plan_dft_3d( PLAN_BACKWARD, Nx,Ny,Nz, zdata, &
                            zdata, FFTW_BACKWARD, FFTW_ESTIMATE )
    CALL dfftw_execute_dft( PLAN_BACKWARD, zdata, zdata )
    CALL dfftw_destroy_plan( PLAN_BACKWARD )
  !
  ELSE ! forward transform
    CALL dfftw_plan_dft_3d( PLAN_FORWARD, Nx, Ny, Nz, zdata, &
                            zdata, FFTW_FORWARD,FFTW_ESTIMATE )
    CALL dfftw_execute_dft( PLAN_FORWARD, zdata, zdata )
    ! Scale the result
    zdata = zdata / ( Nx*Ny*Nz ) !TODO: Use zdscal
    CALL dfftw_destroy_plan( PLAN_FORWARD )
  ENDIF

END SUBROUTINE

