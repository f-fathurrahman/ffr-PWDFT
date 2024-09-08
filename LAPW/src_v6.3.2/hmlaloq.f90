
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlaloq(ias,ngp,ngpq,apwalm,apwalmq,ld,hq)
use modmain
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: hq(ld,*)
! local variables
integer is,io,ilo
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3,i,j
complex(8) z1
is=idxis(ias)
do ilo=1,nlorb(is)
  l1=lorbl(ilo,is)
  do m1=-l1,l1
    lm1=idxlm(l1,m1)
    lm3=0
    do l3=0,lmaxapw
      do m3=-l3,l3
        lm3=lm3+1
        do io=1,apword(l3,is)
          z1=0.d0
          do l2=0,lmaxo
            if (mod(l1+l2+l3,2).eq.0) then
              do m2=-l2,l2
                lm2=idxlm(l2,m2)
                z1=z1+gntyry(lm2,lm3,lm1)*hloa(lm2,io,l3,ilo,ias)
              end do
            end if
          end do
          if (abs(dble(z1))+abs(aimag(z1)).gt.1.d-14) then
            j=ngp+idxlo(lm1,ilo,ias)
            do i=1,ngpq
              hq(i,j)=hq(i,j)+conjg(z1*apwalmq(i,io,lm3))
            end do
            i=ngpq+idxlo(lm1,ilo,ias)
            do j=1,ngp
              hq(i,j)=hq(i,j)+z1*apwalm(j,io,lm3)
            end do
          end if
        end do
      end do
    end do
  end do
end do
return
end subroutine

