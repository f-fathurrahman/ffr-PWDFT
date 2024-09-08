
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findprimcell
! !INTERFACE:
subroutine findprimcell
! !USES:
use modmain
! !DESCRIPTION:
!   This routine finds the smallest primitive cell which produces the same
!   crystal structure as the conventional cell. This is done by searching
!   through all the vectors which connect atomic positions and finding those
!   which leave the crystal structure invariant. Of these, the three shortest
!   which produce a non-zero unit cell volume are chosen.
!
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,js,ia,ja,ka,na
integer i1,i2,i3,i,j,n
real(8) v1(3),v2(3),v3(3)
real(8) t1,t2
! allocatable arrays
real(8), allocatable :: dp(:),vp(:,:)
do is=1,nspecies
  do ia=1,natoms(is)
! make sure all atomic coordinates are in [0,1)
    call r3frac(epslat,atposl(:,ia,is))
! determine atomic Cartesian coordinates
    call r3mv(avec,atposl(:,ia,is),atposc(:,ia,is))
  end do
end do
! find the smallest set of atoms
is=1
do js=1,nspecies
! if a species has only one atom the cell must be primitive
  if (natoms(js).eq.1) return
  if (natoms(js).lt.natoms(is)) is=js
end do
n=27*natoms(is)
allocate(dp(n),vp(3,n))
! generate set of possible lattice vectors
n=0
do ia=1,natoms(is)
  v1(:)=atposl(:,ia,is)-atposl(:,1,is)
  do i1=-1,1
    v2(1)=v1(1)+dble(i1)
    do i2=-1,1
      v2(2)=v1(2)+dble(i2)
      do i3=-1,1
        v2(3)=v1(3)+dble(i3)
        t1=abs(v2(1))+abs(v2(2))+abs(v2(3))
        if (t1.lt.epslat) goto 20
! check if vector v2 leaves conventional cell invariant
        do js=1,nspecies
          do ja=1,natoms(js)
            v3(:)=atposl(:,ja,js)+v2(:)
            call r3frac(epslat,v3)
            do ka=1,natoms(js)
! check both positions and magnetic fields are the same
              t1=sum(abs(atposl(:,ka,js)-v3(:)))
              t2=sum(abs(bfcmt0(:,ja,js)-bfcmt0(:,ka,js)))
              if ((t1.lt.epslat).and.(t2.lt.epslat)) goto 10
            end do
! atom ja has no equivalent under translation by v2
            goto 20
10 continue
          end do
        end do
! cell invariant under translation by v2, so add to list
        n=n+1
        call r3mv(avec,v2,vp(:,n))
        dp(n)=sqrt(vp(1,n)**2+vp(2,n)**2+vp(3,n)**2)
20 continue
      end do
    end do
  end do
end do
if (n.eq.0) then
  write(*,*)
  write(*,'("Error(findprimcell): cannot find any lattice vectors")')
  write(*,*)
  stop
end if
! find the shortest lattice vector
j=1
t1=1.d8
do i=1,n
  if (dp(i).lt.t1+epslat) then
    j=i
    t1=dp(i)
  end if
end do
avec(:,1)=vp(:,j)
! find the next shortest lattice vector not parallel to the first
j=1
t1=1.d8
do i=1,n
  call r3cross(avec(:,1),vp(:,i),v1)
  t2=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
  if (t2.gt.epslat) then
    if (dp(i).lt.t1+epslat) then
      j=i
      t1=dp(i)
    end if
  end if
end do
avec(:,2)=vp(:,j)
! find the next shortest lattice vector which gives non-zero unit cell volume
call r3cross(avec(:,1),avec(:,2),v1)
j=1
t1=1.d8
do i=1,n
  t2=dot_product(vp(:,i),v1(:))
  if (abs(t2).gt.epslat) then
    if (dp(i).lt.t1+epslat) then
      j=i
      t1=dp(i)
    end if
  end if
end do
avec(:,3)=vp(:,j)
call r3minv(avec,ainv)
! remove redundant atoms
do is=1,nspecies
  na=0
  do ia=1,natoms(is)
    call r3mv(ainv,atposc(:,ia,is),v1)
    call r3frac(epslat,v1)
    do ja=1,na
      t1=sum(abs(atposl(:,ja,is)-v1(:)))
      if (t1.lt.epslat) goto 30
    end do
    na=na+1
    atposl(:,na,is)=v1(:)
    call r3mv(avec,atposl(:,na,is),atposc(:,na,is))
! re-index external magnetic fields
    bfcmt0(:,na,is)=bfcmt0(:,ia,is)
! re-index fixed spin moment vectors
    mommtfix(:,na,is)=mommtfix(:,ia,is)
30 continue
  end do
  natoms(is)=na
end do
deallocate(dp,vp)
return
end subroutine
!EOC

