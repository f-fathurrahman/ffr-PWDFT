
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpaa(tor,ias,ngp,apwalm,ld,o)
use modmain
implicit none
! arguments
logical, intent(in) :: tor
integer, intent(in) :: ias,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: o(*)
! local variables
integer is,lmo,io
integer l,m,lm,i
! allocatable arrays
complex(8), allocatable :: a(:,:)
is=idxis(ias)
lmo=lmoapw(is)
allocate(a(lmo,ngp))
i=0
lm=0
do l=0,lmaxapw
  do m=-l,l
    lm=lm+1
    do io=1,apword(l,is)
      i=i+1
      a(i,1:ngp)=apwalm(1:ngp,io,lm)
    end do
  end do
end do
call zmctmu(tor,lmo,ngp,a,a,ld,o)
deallocate(a)
return
end subroutine

