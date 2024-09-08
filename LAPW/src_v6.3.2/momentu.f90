
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine momentu
use modmain
use modulr
implicit none
! local variables
integer ifq,idm,is,ias
integer nrc,nrci
real(8) t1
! allocatable arrays
real(8), allocatable :: rfft(:)
complex(8), allocatable :: zfft(:)
! external functions
real(8) ddot
complex(8) zfmtint
external ddot,zfmtint
if (.not.spinpol) return
allocate(rfft(nqpt),zfft(nfqrz))
! calculate muffin-tin moments
mommttot(:)=0.d0
do idm=1,ndmag
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    do ifq=1,nfqrz
      zfft(ifq)=zfmtint(nrc,nrci,wrcmt(:,is),magqmt(:,ias,idm,ifq))
    end do
    mommt(idm,ias)=dble(zfft(1))
    mommttot(idm)=mommttot(idm)+mommt(idm,ias)
    call rzfftifc(3,ngridq,1,rfft,zfft)
    mommtru(idm,ias,:)=rfft(:)
  end do
end do
! calculate interstitial moment
do idm=1,ndmag
  t1=ddot(ngtot,magqir(:,idm,1),2,cfunir,1)
  momir(idm)=t1*omega/dble(ngtot)
end do
momtot(:)=mommttot(:)+momir(:)
! total moment magnitude
if (ncmag) then
  momtotm=sqrt(momtot(1)**2+momtot(2)**2+momtot(3)**2)
else
  momtotm=abs(momtot(1))
end if
deallocate(rfft,zfft)
return
end subroutine

