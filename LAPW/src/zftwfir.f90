
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zftwfir(ngp,igpig,wfir)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(inout) :: wfir(ngtc,nspinor,nstsv)
! local variables
integer ist,ispn,jspn
integer igp,ifg
real(8) t0
! allocatable arrays
complex(8), allocatable :: z(:)

t0=1.d0/sqrt(omega)
allocate(z(ngkmax))
do ist=1,nstsv
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    call zcopy(ngp(jspn),wfir(:,ispn,ist),1,z,1)
    wfir(:,ispn,ist)=0.d0
    do igp=1,ngp(jspn)
      ifg=igfc(igpig(igp,jspn))
      wfir(ifg,ispn,ist)=t0*z(igp)
    end do
    call zfftifc(3,ngdc,1,wfir(:,ispn,ist))
  end do
end do
deallocate(z)
return
end subroutine

