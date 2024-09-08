
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: spline
! !INTERFACE:
subroutine spline(n,x,f,cf)
! !INPUT/OUTPUT PARAMETERS:
!   n  : number of points (in,integer)
!   x  : abscissa array (in,real(n))
!   f  : input data array (in,real(n))
!   cf : cubic spline coefficients (out,real(3,n))
! !DESCRIPTION:
!   Calculates the coefficients of a cubic spline fitted to input data. In other
!   words, given a set of data points $f_i$ defined at $x_i$, where
!   $i=1\ldots n$, the coefficients $c_j^i$ are determined such that
!   $$ y_i(x)=f_i+c_1^i(x-x_i)+c_2^i(x-x_i)^2+c_3^i(x-x_i)^3, $$
!   is the interpolating function for $x\in[x_i,x_{i+1})$. The coefficients are
!   determined piecewise by fitting a cubic polynomial to adjacent points.
!
! !REVISION HISTORY:
!   Created November 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n),f(n)
real(8), intent(out) :: cf(3,n)
! local variables
integer i
real(8) x0,x1,x2,x3,y0,y1,y2,y3
real(8) t0,t1,t2,t3,t4,t5,t6,t7
if (n.le.0) then
  write(*,*)
  write(*,'("Error(spline): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
if (n.eq.1) then
  cf(:,1)=0.d0
  return
end if
if (n.eq.2) then
  cf(1,1)=(f(2)-f(1))/(x(2)-x(1))
  cf(2:3,1)=0.d0
  cf(1,2)=cf(1,1)
  cf(2:3,2)=0.d0
  return
end if
if (n.eq.3) then
  x0=x(1)
  x1=x(2)-x0; x2=x(3)-x0
  y0=f(1)
  y1=f(2)-y0; y2=f(3)-y0
  t0=1.d0/(x1*x2*(x2-x1))
  t3=x1*y2; t4=x2*y1
  t1=t0*(x2*t4-x1*t3)
  t2=t0*(t3-t4)
  cf(1,1)=t1
  cf(2,1)=t2
  cf(3,1)=0.d0
  t3=2.d0*t2
  cf(1,2)=t1+t3*x1
  cf(2,2)=t2
  cf(3,2)=0.d0
  cf(1,3)=t1+t3*x2
  cf(2,3)=t2
  cf(3,3)=0.d0
  return
end if
x0=x(1)
x1=x(2)-x0; x2=x(3)-x0; x3=x(4)-x0
t4=x1-x2; t5=x1-x3; t6=x2-x3
y0=f(1)
y1=f(2)-y0; y2=f(3)-y0; y3=f(4)-y0
t1=x1*x2*y3; t2=x2*x3*y1; t3=x1*x3
t0=1.d0/(x2*t3*t4*t5*t6)
t3=t3*y2
t7=t0*(t1*t4+t2*t6-t3*t5)
t4=x1**2; t5=x2**2; t6=x3**2
y1=t3*t6-t1*t5; y3=t2*t5-t3*t4; y2=t1*t4-t2*t6
t1=t0*(x1*y1+x2*y2+x3*y3)
t2=-t0*(y1+y2+y3)
cf(1,1)=t1; cf(2,1)=t2; cf(3,1)=t7
cf(1,2)=t1+2.d0*t2*x1+3.d0*t7*t4
cf(2,2)=t2+3.d0*t7*x1
cf(3,2)=t7
if (n.eq.4) then
  cf(1,3)=t1+2.d0*t2*x2+3.d0*t7*t5
  cf(2,3)=t2+3.d0*t7*x2
  cf(3,3)=t7
  cf(1,4)=t1+2.d0*t2*x3+3.d0*t7*t6
  cf(2,4)=t2+3.d0*t7*x3
  cf(3,4)=t7
  return
end if
do i=3,n-2
  x0=x(i)
  x1=x(i-1)-x0; x2=x(i+1)-x0; x3=x(i+2)-x0
  t4=x1-x2; t5=x1-x3; t6=x2-x3
  y0=f(i)
  y1=f(i-1)-y0; y2=f(i+1)-y0; y3=f(i+2)-y0
  t1=x1*x2*y3; t2=x2*x3*y1; t3=x1*x3
  t0=1.d0/(x2*t3*t4*t5*t6)
  t3=t3*y2
  t7=t0*(t1*t4+t2*t6-t3*t5)
  t4=x1**2; t5=x2**2; t6=x3**2
  y1=t3*t6-t1*t5; y2=t1*t4-t2*t6; y3=t2*t5-t3*t4
  t1=t0*(x1*y1+x2*y2+x3*y3)
  t2=-t0*(y1+y2+y3)
  cf(1,i)=t1; cf(2,i)=t2; cf(3,i)=t7
end do
cf(1,n-1)=t1+2.d0*t2*x2+3.d0*t7*t5
cf(2,n-1)=t2+3.d0*t7*x2
cf(3,n-1)=t7
cf(1,n)=t1+2.d0*t2*x3+3.d0*t7*t6
cf(2,n)=t2+3.d0*t7*x3
cf(3,n)=t7
return
end subroutine
!EOC

