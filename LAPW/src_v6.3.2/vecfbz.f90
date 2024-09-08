
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: vecfbz
! !INTERFACE:
subroutine vecfbz(eps,bvec,vpl)
! !INPUT/OUTPUT PARAMETERS:
!   eps  : zero component tolerance (in,real)
!   bvec : reciprocal lattice vectors (in,real(3,3))
!   vpl  : input vector in lattice coordinates (inout,real(3))
! !DESCRIPTION:
!   Maps a vector in lattice coordinates to the first Brillouin zone. This is
!   done by first removing its integer components and then adding primitive
!   reciprocal lattice vectors until the shortest vector is found.
!
! !REVISION HISTORY:
!   Created September 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps,bvec(3,3)
real(8), intent(inout) :: vpl(3)
! local variables
integer i1,i2,i3,j1,j2,j3
real(8) v0(3),v1(3),v2(3),v3(3),t1,t2
! map vector to [0,1) interval
call r3frac(eps,vpl)
v0(:)=bvec(:,1)*vpl(1)+bvec(:,2)*vpl(2)+bvec(:,3)*vpl(3)
t1=v0(1)**2+v0(2)**2+v0(3)**2
j1=0; j2=0; j3=0
do i1=-1,0
  v1(:)=v0(:)+dble(i1)*bvec(:,1)
  do i2=-1,0
    v2(:)=v1(:)+dble(i2)*bvec(:,2)
    do i3=-1,0
      v3(:)=v2(:)+dble(i3)*bvec(:,3)
      t2=v3(1)**2+v3(2)**2+v3(3)**2
      if (t2.lt.t1+eps) then
        j1=i1; j2=i2; j3=i3
        t1=t2
      end if
    end do
  end do
end do
vpl(1)=vpl(1)+dble(j1)
vpl(2)=vpl(2)+dble(j2)
vpl(3)=vpl(3)+dble(j3)
return
end subroutine
!EOC

