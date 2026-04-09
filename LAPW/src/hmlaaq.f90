
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlaaq(ias,ngp,ngpq,apwalm,apwalmq,ld,hq)
use modmain
implicit none
integer, intent(in) :: ias
integer, intent(in) :: ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: hq(*)
! local variables
integer is,lmo,io,jo,i
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3
real(8) t0
complex(8) z1
! allocatable arrays
complex(8), allocatable :: a(:,:),b(:,:)
is=idxis(ias)
lmo=lmoapw(is)
allocate(a(lmo,ngpq),b(lmo,ngp))
t0=0.5d0*rmt(is)**2
i=0
lm1=0
do l1=0,lmaxapw
  do m1=-l1,l1
    lm1=lm1+1
    do io=1,apword(l1,is)
      i=i+1
      b(i,:)=0.d0
      lm3=0
      do l3=0,lmaxapw
        do m3=-l3,l3
          lm3=lm3+1
          do jo=1,apword(l3,is)
            z1=0.d0
            do l2=0,lmaxo
              if (mod(l1+l2+l3,2).eq.0) then
                do m2=-l2,l2
                  lm2=idxlm(l2,m2)
                  z1=z1+gntyry(lm2,lm3,lm1)*haa(lm2,jo,l3,io,l1,ias)
                end do
              end if
            end do
            if (abs(dble(z1))+abs(aimag(z1)).gt.1.d-14) then
              call zaxpy(ngp,z1,apwalm(:,jo,lm3),1,b(i,1),lmo)
            end if
          end do
        end do
      end do
! kinetic surface contribution
      do jo=1,apword(l1,is)
        z1=t0*apwfr(nrmt(is),1,io,l1,ias)*apwdfr(jo,l1,ias)
        call zaxpy(ngp,z1,apwalm(:,jo,lm1),1,b(i,1),lmo)
      end do
      a(i,1:ngpq)=apwalmq(1:ngpq,io,lm1)
    end do
  end do
end do
call zmctm(lmo,ngpq,ngp,a,b,ld,hq)
deallocate(a,b)
return
end subroutine

