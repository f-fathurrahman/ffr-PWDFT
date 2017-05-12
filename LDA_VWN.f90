!! PURPOSE:
!!
!!   These subroutines calculate XC energy per volume (`SUBROUTINE excVWN`)
!!   and its derivative w.r.t density (`SUBROUTINE excpVWN`).
!!   for VWN parameterization of the exchange correlation energy.
!!   Based on T.A. Arias code for his DFT course.
!!
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! NOTES:
!! 
!!   This subroutine should be used for testing purpose only.
!!   In the future, its functionality should be replaced by call to LibXC.

SUBROUTINE excVWN( Npts, rho, epsxc )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  INTEGER :: Npts
  REAL(8) :: rho(Npts)
  REAL(8) :: epsxc(Npts)
  !
  REAL(8), ALLOCATABLE :: rs(:), x(:), XX(:)
  ! Constants
  REAL(8), PARAMETER :: X1 = 0.75*(3.0/(2.0*PI))**(2.0/3.0)
  REAL(8), PARAMETER :: A = 0.0310907
  REAL(8), PARAMETER :: x0 = -0.10498
  REAL(8), PARAMETER :: b = 3.72744
  REAL(8), PARAMETER :: c = 12.9352
  REAL(8) :: Q, XX0
  !
  ALLOCATE( rs(Npts), x(Npts), XX(Npts) )
  !
  Q = sqrt(4.0*c - b*b)
  XX0 = x0*x0 + b*x0 + c
  !
  rs = (4.0*PI/3*rho)**(-1.0/3.0) ! Added internal conversion to rs
  x = sqrt(rs)
  XX = x*x + b*x + c
  !
  epsxc = -X1/rs + A*( log(x*x/XX) + 2*b/Q*atan(Q/(2*x + b) ) &
        -(b*x0)/XX0*( log( (x-x0)*(x-x0)/XX ) + 2.0*(2.0*x0 + b)/Q*atan(Q/(2*x + b)) ))
  !
  DEALLOCATE( rs, x, XX )

END SUBROUTINE 


SUBROUTINE excpVWN( Npts, rho, depsxc )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  INTEGER :: Npts
  REAL(8) :: rho(Npts)
  REAL(8) :: depsxc(Npts)
  ! Constants
  REAL(8), PARAMETER :: X1 = 0.75d0*(3.d0/(2.d0*PI))**(2.d0/3.d0)
  REAL(8), PARAMETER :: A = 0.0310907d0
  REAL(8), PARAMETER :: x0 = -0.10498d0
  REAL(8), PARAMETER :: b = 3.72744d0 
  REAL(8), PARAMETER :: c = 12.9352
  REAL(8) :: Q, XX0
  REAL(8), ALLOCATABLE :: rs(:), x(:), XX(:), dx(:)

  ALLOCATE( rs(Npts), x(Npts), XX(Npts), dx(Npts) )

  Q = sqrt(4.*c - b*b)
  XX0 = x0*x0 + b*x0 + c
  rs = (4.d0*pi/3.d0*rho)**(-1.d0/3.d0) ! Added internal conversion to rs
  x = sqrt(rs)
  XX = x*x + b*x + c
  dx = 0.5d0/x ! Chain rule needs dx/drho!
  depsxc = dx*( 2.d0*X1 / (rs*x) &
        + A*( 2d0/x - (2d0*x + b)/XX- 4d0*b/(Q*Q + (2d0*x + b) * (2d0*x + b)) &
        - (b*x0)/XX0*( 2.d0/(x-x0) - (2d0*x + b)/XX - 4d0*(2d0*x0 + b) / &
          ( Q*Q + (2d0*x + b) * (2d0*x + b)) ) ))
  depsxc = (-rs/(3d0*rho)) * depsxc ! Added d(rs)/dn from chain rule from rs to n conv

  DEALLOCATE( rs, x, XX, dx )
END

