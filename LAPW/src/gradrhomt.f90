
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradrhomt
use modmain
use modphonon
implicit none
! local variables
integer nr,nri,np
! allocatable arrays
complex(8), allocatable :: zfmt(:),gzfmt(:,:)
! add gradient contribution from rigid shift of muffin-tin
nr=nrmt(isph)
nri=nrmti(isph)
np=npmt(isph)
allocate(zfmt(np),gzfmt(np,3))
! convert the density to complex spherical harmonic expansion
call rtozfmt(nr,nri,rhomt(:,iasph),zfmt)
! compute the gradient
call gradzfmt(nr,nri,rlmt(:,1,isph),rlmt(:,-1,isph),zfmt,np,gzfmt)
! subtract from the density derivative
drhomt(1:np,iasph)=drhomt(1:np,iasph)-gzfmt(1:np,ipph)
deallocate(zfmt,gzfmt)
return
end subroutine

