
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine chargeu
use modmain
use modulr
implicit none
! local variables
integer ifq,is,ias
integer nrc,nrci
real(8) t1
! allocatable arrays
real(8), allocatable :: rfft(:)
complex(8), allocatable :: zfft(:)
! external functions
real(8) ddot
complex(8) zfmtint
external ddot,zfmtint
allocate(rfft(nqpt),zfft(nfqrz))
! calculate muffin-tin charges
chgmttot=0.d0
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  do ifq=1,nfqrz
    zfft(ifq)=zfmtint(nrc,nrci,wrcmt(:,is),rhoqmt(:,ias,ifq))
  end do
  chgmt(ias)=dble(zfft(1))
  chgmttot=chgmttot+chgmt(ias)
  call rzfftifc(3,ngridq,1,rfft,zfft)
  chgmtru(ias,:)=rfft(:)
end do
! calculate interstitial charge
t1=ddot(ngtot,rhoqir(:,1),2,cfunir,1)
chgir=t1*omega/dble(ngtot)
! total calculated charge
chgcalc=chgmttot+chgir
deallocate(rfft,zfft)
return
end subroutine

