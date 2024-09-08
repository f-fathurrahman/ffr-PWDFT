
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfirsm(m,rfir)
use modmain
implicit none
! arguments
integer, intent(in) :: m
real(8), intent(inout) :: rfir(ngtot)
! local variables
integer ig,ifg
real(8) t0,t1,t2
! allocatable arrays
complex(8), allocatable :: zfft(:)
if (m.le.0) return
allocate(zfft(ngtot))
zfft(:)=rfir(:)
call zfftifc(3,ngridg,-1,zfft)
t0=dble(2*m)
t1=1.d0/gmaxvr
do ig=1,ngtot
  ifg=igfft(ig)
  t2=t1*gc(ig)
  t2=exp(-t0*t2**4)
  zfft(ifg)=t2*zfft(ifg)
end do
call zfftifc(3,ngridg,1,zfft)
rfir(:)=dble(zfft(:))
deallocate(zfft)
return
end subroutine

