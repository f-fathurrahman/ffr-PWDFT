
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writestulr
use modmain
use modulr
implicit none
! local variables
integer ifq,idm
open(40,file='STATE_ULR.OUT',form='UNFORMATTED')
write(40) version
write(40) natmtot
write(40) npcmtmax
write(40) ngtot
write(40) ndmag
write(40) fsmtype
write(40) nqpt
write(40) nfqrz
write(40) ivq
write(40) iqrzf
! write the ultra long-range density in Q-space
do ifq=1,nfqrz
  write(40) rhoqmt(:,:,ifq)
  write(40) rhoqir(:,ifq)
end do
! write the Kohn-Sham effective potential in Q-space
do ifq=1,nfqrz
  write(40) vsqmt(:,:,ifq)
  write(40) vsqir(:,ifq)
end do
! write the external Coulomb potential in Q-space
do ifq=1,nfqrz
  write(40) vclq(ifq)
end do
if (spinpol) then
! write the magnetisation in Q-space
  do ifq=1,nfqrz
    do idm=1,ndmag
      write(40) magqmt(:,:,idm,ifq)
      write(40) magqir(:,idm,ifq)
    end do
  end do
! write the Kohn-Sham effective magnetic field in Q-space
  do ifq=1,nfqrz
    do idm=1,ndmag
      write(40) bsqmt(:,:,idm,ifq)
      write(40) bsqir(:,idm,ifq)
    end do
  end do
! write the external magnetic fields in Q-space
  do ifq=1,nfqrz
    do idm=1,ndmag
      write(40) bfcq(idm,ifq)
      write(40) bfcmtq(:,idm,ifq)
    end do
  end do
! write fixed spin moment magnetic fields
  if (fsmtype.ne.0) then
    write(40) bfsmc
    write(40) bfsmcmt
  end if
end if
close(40)
return
end subroutine

