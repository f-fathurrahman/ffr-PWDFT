
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readstulr
use modmain
use modulr
implicit none
! local variables
integer iq,jq,ifq,jfq
integer idm,i1,i2,i3
integer version_(3),ios
integer natmtot_,npcmtmax_,ngtot_
integer ndmag_,fsmtype_,nqpt_,nfqrz_
complex(8) z1
! automatic arrays
complex(8) zv(natmtot)
! allocatable arrays
integer, allocatable :: ivq_(:,:),iqrzf_(:),map(:)
complex(8), allocatable :: zfmt(:,:),zfir(:)
open(40,file='STATE_ULR.OUT',form='UNFORMATTED',status='OLD',iostat=ios)
if (ios.ne.0) then
  write(*,*)
  write(*,'("Error(readstulr): error opening STATE_ULR.OUT")')
  write(*,*)
  stop
end if
read(40) version_
if ((version(1).ne.version_(1)).or.(version(2).ne.version_(2)).or. &
    (version(3).ne.version_(3))) then
  write(*,*)
  write(*,'("Warning(readstulr): different versions")')
  write(*,'(" current       : ",I3.3,".",I3.3,".",I3.3)') version
  write(*,'(" STATE_ULR.OUT : ",I3.3,".",I3.3,".",I3.3)') version_
end if
read(40) natmtot_
if (natmtot.ne.natmtot_) then
  write(*,*)
  write(*,'("Error(readstulr): differing natmtot")')
  write(*,'(" current       : ",I6)') natmtot
  write(*,'(" STATE_ULR.OUT : ",I6)') natmtot_
  write(*,*)
  stop
end if
read(40) npcmtmax_
if (npcmtmax.ne.npcmtmax_) then
  write(*,*)
  write(*,'("Error(readstulr): differing npcmtmax")')
  write(*,'(" current       : ",I6)') npcmtmax
  write(*,'(" STATE_ULR.OUT : ",I6)') npcmtmax_
  write(*,*)
  stop
end if
read(40) ngtot_
if (ngtot.ne.ngtot_) then
  write(*,*)
  write(*,'("Error(readstulr): differing ngtot")')
  write(*,'(" current       : ",I8)') ngtot
  write(*,'(" STATE_ULR.OUT : ",I8)') ngtot_
  write(*,*)
  stop
end if
read(40) ndmag_
if (ndmag.ne.ndmag_) then
  write(*,*)
  write(*,'("Error(readstulr): differing ndmag")')
  write(*,'(" current       : ",I1)') ndmag
  write(*,'(" STATE_ULR.OUT : ",I1)') ndmag_
  write(*,*)
  stop
end if
read(40) fsmtype_
if (fsmtype.ne.fsmtype_) then
  write(*,*)
  write(*,'("Error(readstulr): differing fsmtype")')
  write(*,'(" current       : ",I4)') fsmtype
  write(*,'(" STATE_ULR.OUT : ",I4)') fsmtype_
  write(*,*)
  stop
end if
read(40) nqpt_
if (nqpt_.le.0) then
  write(*,*)
  write(*,'("Error(readstulr): nqpt_ <= 0 : ",I8)') nqpt_
  write(*,*)
  stop
end if
read(40) nfqrz_
if (nfqrz_.le.0) then
  write(*,*)
  write(*,'("Error(readstulr): nfqrz_ <= 0 : ",I8)') nfqrz_
  write(*,*)
  stop
end if
allocate(ivq_(3,nqpt_),iqrzf_(nfqrz_),map(nfqrz_))
read(40) ivq_
read(40) iqrzf_
! generate map from old Q-vector grid to new
map(:)=0
do ifq=1,nfqrz_
  iq=iqrzf_(ifq)
  i1=ivq_(1,iq); i2=ivq_(2,iq); i3=ivq_(3,iq)
  if ((i1.ge.intq(1,1)).and.(i1.le.intq(2,1)).and. &
      (i2.ge.intq(1,2)).and.(i2.le.intq(2,2)).and. &
      (i3.ge.intq(1,3)).and.(i3.le.intq(2,3))) then
    jq=ivqiq(i1,i2,i3)
    jfq=ifqrz(jq)
    map(ifq)=jfq
  end if
end do
deallocate(ivq_,iqrzf_)
allocate(zfmt(npcmtmax,natmtot),zfir(ngtot))
! read the Q-space density
rhoqmt(:,:,:)=0.d0
rhoqir(:,:)=0.d0
do ifq=1,nfqrz_
  jfq=map(ifq)
  if (jfq.gt.0) then
    read(40) rhoqmt(:,:,jfq)
    read(40) rhoqir(:,jfq)
  else
    read(40) zfmt
    read(40) zfir
  end if
end do
! read the Q-space Kohn-Sham potential
vsqmt(:,:,:)=0.d0
vsqir(:,:)=0.d0
do ifq=1,nfqrz_
  jfq=map(ifq)
  if (jfq.gt.0) then
    read(40) vsqmt(:,:,jfq)
    read(40) vsqir(:,jfq)
  else
    read(40) zfmt
    read(40) zfir
  end if
end do
! read the external Coulomb potential in Q-space
vclq(:)=0.d0
do ifq=1,nfqrz_
  jfq=map(ifq)
  if (jfq.gt.0) then
    read(40) vclq(jfq)
  else
    read(40) z1
  end if
end do
if (spinpol) then
! read the Q-space magnetisation density
  magqmt(:,:,:,:)=0.d0
  magqir(:,:,:)=0.d0
  do ifq=1,nfqrz_
    jfq=map(ifq)
    if (jfq.gt.0) then
      do idm=1,ndmag
        read(40) magqmt(:,:,idm,jfq)
        read(40) magqir(:,idm,jfq)
      end do
    else
      do idm=1,ndmag
        read(40) zfmt
        read(40) zfir
      end do
    end if
  end do
  bsqmt(:,:,:,:)=0.d0
  bsqir(:,:,:)=0.d0
  do ifq=1,nfqrz_
    jfq=map(ifq)
    if (jfq.gt.0) then
      do idm=1,ndmag
        read(40) bsqmt(:,:,idm,jfq)
        read(40) bsqir(:,idm,jfq)
      end do
    else
      do idm=1,ndmag
        read(40) zfmt
        read(40) zfir
      end do
    end if
  end do
! read the external magnetic fields in Q-space
  bfcq(:,:)=0.d0
  bfcmtq(:,:,:)=0.d0
  do ifq=1,nfqrz_
    jfq=map(ifq)
    if (jfq.gt.0) then
      do idm=1,ndmag
        read(40) bfcq(idm,jfq)
        read(40) bfcmtq(:,idm,jfq)
      end do
    else
      do idm=1,ndmag
        read(40) z1
        read(40) zv
      end do
    end if
  end do
! read fixed spin moment effective fields
  if (fsmtype.ne.0) then
    read(40) bfsmc
    read(40) bfsmcmt
  end if
end if
close(40)
deallocate(map,zfmt,zfir)
return
end subroutine

