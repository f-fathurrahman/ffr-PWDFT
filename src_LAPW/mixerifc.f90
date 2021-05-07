
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mixerifc(mtype,n,v,dv,nwork,work)
use modmain
implicit none
! arguments
integer, intent(in) :: mtype,n
real(8), intent(inout) :: v(n)
real(8), intent(out) :: dv
integer, intent(inout) :: nwork
real(8), intent(inout) :: work(*)
! local variables
select case(mtype)
case(0)
! linear mixing
  if (nwork.le.0) then
    nwork=n
    return
  end if
  call mixlinear(iscl,amixpm(1),n,v,work,dv)
case(1)
! adaptive linear mixing
  if (nwork.le.0) then
    nwork=3*n
    return
  end if
  call mixadapt(iscl,amixpm(1),amixpm(2),n,v,work,work(n+1),work(2*n+1),dv)
case(3)
! Broyden mixing
  if (nwork.le.0) then
    nwork=(4+2*mixsdb)*n+mixsdb**2
    return
  end if
  call mixbroyden(iscl,n,mixsdb,broydpm(1),broydpm(2),v,work,work(2*n+1), &
   work(4*n+1),work((4+mixsdb)*n+1),work((4+2*mixsdb)*n+1),dv)
case default
  write(*,*)
  write(*,'("Error(mixerifc): mtype not defined : ",I8)') mtype
  write(*,*)
  stop
end select
return
end subroutine

subroutine getmixdata(mtype,mixdescr)
implicit none
! arguments
integer, intent(in) :: mtype
character(*), intent(out) :: mixdescr
select case(mtype)
case(0)
  mixdescr='Linear mixing'
case(1)
  mixdescr='Adaptive linear mixing'
case(3)
  mixdescr='Broyden mixing, J. Phys. A: Math. Gen. 17, L317 (1984)'
case default
  write(*,*)
  write(*,'("Error(getmixdata): mixtype not defined : ",I8)') mtype
  write(*,*)
  stop
end select
return
end subroutine

