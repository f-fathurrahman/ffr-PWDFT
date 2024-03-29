
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rotaxang
! !INTERFACE:
subroutine rotaxang(eps,rot,det,v,th)
! !INPUT/OUTPUT PARAMETERS:
!   eps : zero vector tolerance (in,real)
!   rot : rotation matrix (in,real(3,3))
!   det : matrix determinant (out,real)
!   v   : normalised axis vector (out,real(3))
!   th  : rotation angle (out,real)
! !DESCRIPTION:
!   Given a rotation matrix
!   $$ R(\hat{\bf v},\theta)=
!     \left(\begin{matrix}
!     \cos\theta+x^2(1-\cos\theta) &
!     xy(1-\cos\theta)+z\sin\theta &
!     xz(1-\cos\theta)-y\sin\theta \\
!     xy(1-\cos\theta)-z\sin\theta &
!     \cos\theta+y^2(1-\cos\theta) &
!     yz(1-\cos\theta)+x\sin\theta \\
!     xz(1-\cos\theta)+y\sin\theta &
!     yz(1-\cos\theta)-x\sin\theta &
!     \cos\theta+z^2(1-\cos\theta)
!    \end{matrix}\right), $$
!   this routine determines the axis of rotation $\hat{\bf v}$ and the angle of
!   rotation $\theta$. If $R$ corresponds to an improper rotation then only the
!   proper part is used and {\tt det} is set to $-1$. The rotation convention
!   follows the `right-hand rule'.
!
! !REVISION HISTORY:
!   Created December 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(in) :: rot(3,3)
real(8), intent(out) :: det
real(8), intent(out) :: v(3),th
! local variables
real(8), parameter :: pi=3.1415926535897932385d0
real(8) rotp(3,3),t1,t2
! find the determinant
det=rot(1,2)*rot(2,3)*rot(3,1)-rot(1,3)*rot(2,2)*rot(3,1) &
   +rot(1,3)*rot(2,1)*rot(3,2)-rot(1,1)*rot(2,3)*rot(3,2) &
   +rot(1,1)*rot(2,2)*rot(3,3)-rot(1,2)*rot(2,1)*rot(3,3)
if (abs(det-1.d0).lt.eps) then
  det=1.d0
else if (abs(det+1.d0).lt.eps) then
  det=-1.d0
else
  goto 10
end if
! proper rotation matrix
rotp(:,:)=det*rot(:,:)
v(1)=(rotp(2,3)-rotp(3,2))/2.d0
v(2)=(rotp(3,1)-rotp(1,3))/2.d0
v(3)=(rotp(1,2)-rotp(2,1))/2.d0
t1=sqrt(v(1)**2+v(2)**2+v(3)**2)
t2=(rotp(1,1)+rotp(2,2)+rotp(3,3)-1.d0)/2.d0
if (abs(abs(t2)-1.d0).gt.eps) then
! theta not equal to 0 or pi
  th=-atan2(t1,t2)
  v(:)=v(:)/t1
else
! special case of sin(th)=0
  if (t2.gt.0.d0) then
! zero angle: axis arbitrary
    th=0.d0
    v(:)=1.d0/sqrt(3.d0)
  else
! rotation by pi
    th=pi
    if ((rotp(1,1).ge.rotp(2,2)).and.(rotp(1,1).ge.rotp(3,3))) then
      if (rotp(1,1).lt.(-1.d0+eps)) goto 10
      v(1)=sqrt(abs(rotp(1,1)+1.d0)/2.d0)
      v(2)=(rotp(2,1)+rotp(1,2))/(4.d0*v(1))
      v(3)=(rotp(3,1)+rotp(1,3))/(4.d0*v(1))
    else if ((rotp(2,2).ge.rotp(1,1)).and.(rotp(2,2).ge.rotp(3,3))) then
      if (rotp(2,2).lt.(-1.d0+eps)) goto 10
      v(2)=sqrt(abs(rotp(2,2)+1.d0)/2.d0)
      v(3)=(rotp(3,2)+rotp(2,3))/(4.d0*v(2))
      v(1)=(rotp(1,2)+rotp(2,1))/(4.d0*v(2))
    else
      if (rotp(3,3).lt.(-1.d0+eps)) goto 10
      v(3)=sqrt(abs(rotp(3,3)+1.d0)/2.d0)
      v(1)=(rotp(1,3)+rotp(3,1))/(4.d0*v(3))
      v(2)=(rotp(2,3)+rotp(3,2))/(4.d0*v(3))
    end if
  end if
end if
return
10 continue
write(*,*)
write(*,'("Error(rotaxang): invalid rotation matrix:")')
write(*,'(3G18.10)') rot(1,:)
write(*,'(3G18.10)') rot(2,:)
write(*,'(3G18.10)') rot(3,:)
write(*,*)
stop
end subroutine
!EOC

