!------------------------------------------------------------------
SUBROUTINE my_rdiracint(sol, kpa, e, nr, r, vr, nn, g0, g1, f0, f1)
!------------------------------------------------------------------
  IMPLICIT NONE 
  ! arguments
  real(8), intent(in) :: sol
  integer, intent(in) :: kpa
  real(8), intent(in) :: e
  integer, intent(in) :: nr
  real(8), intent(in) :: r(nr),vr(nr)
  integer, intent(out) :: nn
  real(8), intent(out) :: g0(nr),g1(nr)
  real(8), intent(out) :: f0(nr),f1(nr)
  ! local variables
  integer :: ir,ir0
  ! rescaling limit
  real(8), parameter :: rsc=1.d100
  real(8) :: ci,e0,t1,t2,t3,t4
  
  if (nr.lt.4) then
    write(*,*)
    write(*,'("Error(rdiracint): nr < 4 : ",I8)') nr
    write(*,*)
    stop
  end if
  
  ! inverse speed of light
  ci = 1.d0/sol
  
  ! electron rest energy
  e0 = sol**2
  t1 = 2.d0*e0+e
  
  ! determine the r -> 0 boundary values of F and G
  t2 = dble(kpa)/r(1)
  t3 = ci*(t1-vr(1))
  t4 = ci*(vr(1)-e)
  f0(1) = 1.d0
  f1(1) = 0.d0
  g0(1) = (f1(1)-t2*f0(1))/t4
  g1(1) = t3*f0(1)-t2*g0(1)

  ! extrapolate to the first four points
  g1(2:4) = g1(1)
  f1(2:4) = f1(1)
  nn = 0
  do ir = 2,nr
    t2 = dble(kpa)/r(ir)
    t3 = ci*(t1-vr(ir))
    t4 = ci*(vr(ir)-e)
    ir0 = ir-3
    if (ir0 < 1) ir0=1
    g1(ir) = poly3(r(ir0),g1(ir0),r(ir))
    f1(ir) = poly3(r(ir0),f1(ir0),r(ir))
    ! integrate to find wavefunction
    g0(ir) = poly4i(r(ir0),g1(ir0),r(ir))+g0(ir0)
    f0(ir) = poly4i(r(ir0),f1(ir0),r(ir))+f0(ir0)
    ! compute the derivatives
    g1(ir) = t3*f0(ir)-t2*g0(ir)
    f1(ir) = t4*g0(ir)+t2*f0(ir)
    ! integrate for correction
    g0(ir) = poly4i(r(ir0),g1(ir0),r(ir))+g0(ir0)
    f0(ir) = poly4i(r(ir0),f1(ir0),r(ir))+f0(ir0)
    ! compute the derivatives again
    g1(ir) = t3*f0(ir)-t2*g0(ir)
    f1(ir) = t4*g0(ir)+t2*f0(ir)
    ! check for overflow
    if ((abs(g0(ir)).gt.rsc).or.(abs(g1(ir)).gt.rsc).or. &
        (abs(f0(ir)).gt.rsc).or.(abs(f1(ir)).gt.rsc)) then
      ! set the remaining points and return
      g0(ir:nr)=g0(ir)
      g1(ir:nr)=g1(ir)
      f0(ir:nr)=f0(ir)
      f1(ir:nr)=f1(ir)
      return
    end if
    ! check for node
    if (g0(ir-1)*g0(ir).lt.0.d0) nn=nn+1
  end do
return

!--------
contains
!--------

! Local functions: poly3 and poly4i

!------------------------------
real(8) function poly3(xa,ya,x)
!------------------------------
  implicit none
  ! arguments
  real(8) xa(3),ya(3),x
  ! local variables
  real(8) x0,x1,x2,y0,y1,y2
  real(8) c1,c2,t0,t1,t2
  ! evaluate the polynomial coefficients
  x0=xa(1)
  x1=xa(2)-x0; x2=xa(3)-x0
  y0=ya(1)
  y1=ya(2)-y0; y2=ya(3)-y0
  t0=1.d0/(x1*x2*(x2-x1))
  t1=x1*y2; t2=x2*y1
  c1=x2*t2-x1*t1
  c2=t1-t2
  t1=x-x0
  ! evaluate the polynomial
  poly3=y0+t0*t1*(c1+c2*t1)
  return
end function

!-------------------------------
real(8) function poly4i(xa,ya,x)
!-------------------------------
  implicit none
  ! arguments
  real(8), intent(in) :: xa(4),ya(4),x
  ! local variables
  real(8) x0,x1,x2,x3,y0,y1,y2,y3
  real(8) c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
  ! evaluate the polynomial coefficients
  x0=xa(1)
  x1=xa(2)-x0; x2=xa(3)-x0; x3=xa(4)-x0
  y0=ya(1)
  y1=ya(2)-y0; y2=ya(3)-y0; y3=ya(4)-y0
  t4=x1-x2; t5=x1-x3; t6=x2-x3
  t1=x1*x2*y3; t2=x2*x3*y1; t3=x1*x3
  t0=1.d0/(x2*t3*t4*t5*t6)
  t3=t3*y2
  c3=t1*t4+t2*t6-t3*t5
  t4=x1**2; t5=x2**2; t6=x3**2
  c2=t1*(t5-t4)+t2*(t6-t5)+t3*(t4-t6)
  c1=t1*(x2*t4-x1*t5)+t2*(x3*t5-x2*t6)+t3*(x1*t6-x3*t4)
  t1=x-x0
  ! integrate the polynomial
  poly4i=t1*(y0+t0*t1*(0.5d0*c1+t1*(0.3333333333333333333d0*c2+0.25d0*c3*t1)))
  return
end function

end subroutine
