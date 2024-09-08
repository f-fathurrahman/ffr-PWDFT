
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine splintwp(n,wp,f,g)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wp(4,n),f(n)
real(8), intent(out) :: g(n)
! local variables
integer i
g(1)=0.d0
g(2)=wp(1,2)*f(1)+wp(2,2)*f(2)+wp(3,2)*f(3)+wp(4,2)*f(4)
do i=3,n-1
  g(i)=g(i-1)+wp(1,i)*f(i-2)+wp(2,i)*f(i-1)+wp(3,i)*f(i)+wp(4,i)*f(i+1)
end do
g(n)=g(n-1)+wp(1,n)*f(n-3)+wp(2,n)*f(n-2)+wp(3,n)*f(n-1)+wp(4,n)*f(n)
return
end subroutine

