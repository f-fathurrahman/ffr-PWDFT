
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dolpaa(ias,ngp,ngpq,apwalm,apwalmq,dapwalm,dapwalmq,ld,od)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: od(*)
! local variables
integer is,lmo,io
integer l,m,lm,i
! allocatable arrays
complex(8), allocatable :: a(:,:),b(:,:)
if (ias.ne.iasph) return
is=idxis(ias)
lmo=lmoapw(is)
allocate(a(lmo,ngpq),b(lmo,ngp))
i=0
lm=0
do l=0,lmaxapw
  do m=-l,l
    lm=lm+1
    do io=1,apword(l,is)
      i=i+1
      a(i,1:ngpq)=apwalmq(1:ngpq,io,lm)
      b(i,1:ngp)=dapwalm(1:ngp,io,lm)
    end do
  end do
end do
call zmctm(lmo,ngpq,ngp,a,b,ld,od)
i=0
lm=0
do l=0,lmaxapw
  do m=-l,l
    lm=lm+1
    do io=1,apword(l,is)
      i=i+1
      a(i,1:ngpq)=dapwalmq(1:ngpq,io,lm)
      b(i,1:ngp)=apwalm(1:ngp,io,lm)
    end do
  end do
end do
call zmctm(lmo,ngpq,ngp,a,b,ld,od)
deallocate(a,b)
return
end subroutine

