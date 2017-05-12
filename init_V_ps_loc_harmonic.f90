!!
!! Initialize V_ps_loc with harmonic potential.
!! `omega` is the harmonic parameter and `center` is
!! the center of the potential
!!
!! author: Fadjar Fathurrahman
!!
SUBROUTINE init_V_ps_loc_harmonic( omega, center )
  USE m_realspace, ONLY : Npoints, rgrid
  USE m_hamiltonian, ONLY : V_ps_loc
  IMPLICIT NONE 
  INTEGER :: ip
  REAL(8) :: omega
  REAL(8) :: center(3)
  REAL(8) :: dx, dy, dz

  WRITE(*,*)
  WRITE(*,*) 'Initializing V_ps_loc with harmonic potential'
  WRITE(*,*) 'omega = ', omega
  WRITE(*,*) 'center = ', center

  DO ip = 1, Npoints
    dx = rgrid(1,ip) - center(1)
    dy = rgrid(2,ip) - center(2)
    dz = rgrid(3,ip) - center(3)
    V_ps_loc(ip) = 0.5d0 * omega**2 * ( dx**2 + dy**2 + dz**2 )
  ENDDO 

END SUBROUTINE 

