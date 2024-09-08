
! Copyright (C) 2019 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine eveqnbf(n,ld,a,b,w)
implicit none
! arguments
integer, intent(in) :: n,ld
complex(8), intent(inout) :: a(ld,n),b(ld,n)
real(8), intent(out) :: w(n)
! local variables
integer n2,i,j
! allocatable arrays
integer, allocatable :: idx(:)
real(8), allocatable :: w2(:),x(:)
complex(8), allocatable :: h(:,:)
! external functions
real(8) dznrm2
external dznrm2
n2=2*n
! setup the Bogoliubov Hamiltonian
allocate(w2(n2),h(n2,n2))
do j=1,n
  do i=1,n
    h(i,j)=a(i,j)
    h(n+i,n+j)=-conjg(a(i,j))
    h(i,n+j)=b(i,j)
  end do
end do
! find the eigenvalues and eigenvectors
call eveqnz(n2,n2,h,w2)
! select those eigenvectors which have largest U-norms
allocate(idx(n2),x(n2))
do j=1,n2
  x(j)=dznrm2(n,h(:,j),1)
end do
call sortidx(n2,x,idx)
do i=1,n
  j=idx(n2-i+1)
  w(i)=w2(j)
  call zcopy(n,h(1,j),1,a(1,i),1)
  call zcopy(n,h(n+1,j),1,b(1,i),1)
end do
deallocate(idx,w2,x,h)
return
end subroutine

