
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dpotcoul
use modmain
use modphonon
implicit none
! local variables
integer nr,nri,ir,np,i
! allocatable arrays
complex(8), allocatable :: zfmt(:),gzfmt(:,:)
nr=nrmt(isph)
nri=nrmti(isph)
np=npmt(isph)
! solve the complex Poisson's equation in the muffin-tins
call genzvclmt(nrmt,nrmti,nrmtmax,rlmt,wprmt,npmtmax,drhomt,dvclmt)
! calculate the gradient of the monopole potential
allocate(zfmt(npmtmax),gzfmt(npmtmax,3))
zfmt(1:np)=0.d0
i=1
do ir=1,nri
  zfmt(i)=vcln(ir,isph)
  i=i+lmmaxi
end do
do ir=nri+1,nr
  zfmt(i)=vcln(ir,isph)
  i=i+lmmaxo
end do
call gradzfmt(nr,nri,rlmt(:,1,isph),rlmt(:,-1,isph),zfmt,npmtmax,gzfmt)
! subtract gradient component corresponding to the phonon polarisation
dvclmt(1:np,iasph)=dvclmt(1:np,iasph)-gzfmt(1:np,ipph)
deallocate(zfmt,gzfmt)
! solve Poisson's equation in the entire unit cell
call zpotcoul(nrmt,nrmti,npmt,npmti,nrmtmax,rlmt,ngridg,igfft,ngvec,gqc,gclgq, &
 ngvec,jlgqrmt,ylmgq,sfacgq,drhoir,npmtmax,dvclmt,dvclir)
return
end subroutine

