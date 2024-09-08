
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potuinit
use modmain
use modulr
use modrandom
use modomp
implicit none
! local variables
integer ifq,idm,is,ias
real(8) cb,t1
! zero the Kohn-Sham potential for all Q
do ifq=1,nfqrz
  do ias=1,natmtot
    is=idxis(ias)
    vsqmt(1:npcmt(is),ias,ifq)=0.d0
  end do
end do
vsqir(:,:)=0.d0
if (.not.spinpol) return
! zero the Kohn-Sham magnetic field for all Q
do ifq=1,nfqrz
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      bsqmt(1:npcmt(is),ias,idm,ifq)=0.d0
    end do
  end do
end do
bsqir(:,:,:)=0.d0
! coupling constant of the external field (g_e/4c)
cb=gfacte/(4.d0*solsc)
! initialise the external magnetic fields
t1=cb*rndbfcu
do ifq=1,nfqrz
  do idm=1,ndmag
    bfcq(idm,ifq)=t1*cmplx(randomu()-0.5d0,randomu()-0.5d0,8)
    do ias=1,natmtot
      bfcmtq(ias,idm,ifq)=t1*cmplx(randomu()-0.5d0,randomu()-0.5d0,8)
    end do
  end do
end do
bfcq(:,1)=dble(bfcq(:,1))
bfcmtq(:,:,1)=dble(bfcmtq(:,:,1))
return
end subroutine

