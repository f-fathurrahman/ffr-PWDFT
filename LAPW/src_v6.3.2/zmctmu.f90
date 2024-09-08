
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zmctmu(tcr,l,n,a,b,ld,c)
use modomp
implicit none
! arguments
logical, intent(in) :: tcr
integer, intent(in) :: l,n
complex(8), intent(in) :: a(l,n),b(l,n)
integer, intent(in) :: ld
complex(8), intent(out) :: c(*)
! local variables
integer l2,i,j,k,nthd
! external functions
real(8) ddot
complex(8) zdotc
external ddot,zdotc
if (tcr) then
! matrix c is real valued
  l2=2*l
  call holdthd(n,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(k,i) &
!$OMP NUM_THREADS(nthd)
  do j=1,n
    k=(j-1)*ld
    do i=1,j
      k=k+1
      c(k)=c(k)+ddot(l2,a(:,i),1,b(:,j),1)
    end do
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
else
! matrix c is complex valued
  call holdthd(n,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(k,i) &
!$OMP NUM_THREADS(nthd)
  do j=1,n
    k=(j-1)*ld
    do i=1,j
      k=k+1
      c(k)=c(k)+zdotc(l,a(:,i),1,b(:,j),1)
    end do
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
end if
return
end subroutine

