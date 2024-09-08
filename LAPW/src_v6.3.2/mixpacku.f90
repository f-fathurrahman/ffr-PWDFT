
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mixpacku(tpack,n,v)
use modmain
use modulr
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(out) :: n
real(8), intent(inout) :: v(*)
! local variables
integer ifq,idm
n=0
! pack or unpack the Kohn-Sham long-range potential
do ifq=1,nfqrz
  call zfpack(tpack,n,npcmt,npcmtmax,vsqmt(:,:,ifq),vsqir(:,ifq),v)
end do
! pack or unpack the Kohn-Sham magnetic field
if (spinpol) then
  do ifq=1,nfqrz
    do idm=1,ndmag
      call zfpack(tpack,n,npcmt,npcmtmax,bsqmt(:,:,idm,ifq),bsqir(:,idm,ifq),v)
    end do
  end do
end if
return
end subroutine

