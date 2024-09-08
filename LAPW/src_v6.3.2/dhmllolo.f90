
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dhmllolo(ias,ngp,ngpq,ld,dh)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp,ngpq
integer, intent(in) :: ld
complex(8), intent(inout) :: dh(ld,*)
! local variables
integer is,ilo,jlo
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3,i,j
complex(8) z1
is=idxis(ias)
do jlo=1,nlorb(is)
  l3=lorbl(jlo,is)
  do m3=-l3,l3
    lm3=idxlm(l3,m3)
    j=ngp+idxlo(lm3,jlo,ias)
    do ilo=1,nlorb(is)
      l1=lorbl(ilo,is)
      do m1=-l1,l1
        lm1=idxlm(l1,m1)
        i=ngpq+idxlo(lm1,ilo,ias)
        z1=0.d0
        do l2=0,lmaxo
          if (mod(l1+l2+l3,2).eq.0) then
            do m2=-l2,l2
              lm2=idxlm(l2,m2)
              z1=z1+gntyyy(lm2,lm3,lm1)*dhlolo(lm2,jlo,ilo,ias)
            end do
          end if
        end do
        dh(i,j)=dh(i,j)+z1
      end do
    end do
  end do
end do
return
end subroutine

