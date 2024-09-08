
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readdvs(iq,is,ia,ip)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq,is,ia,ip
! local variables
integer js,jas,ios
integer version_(3),nspecies_
integer lmmaxo_,natoms_
integer nrmt_,ngridg_(3)
character(256) fext,fname
! allocatable arrays
complex(8), allocatable :: zfmt(:,:,:)
allocate(zfmt(lmmaxo,nrmtmax,natmtot))
call phfext(iq,is,ia,ip,fext)
fname='DVS'//trim(fext)
open(50,file=trim(fname),form='UNFORMATTED',status='OLD',iostat=ios)
if (ios.ne.0) then
  write(*,*)
  write(*,'("Error(readdvs): error opening ",A)') trim(fname)
  write(*,*)
  stop
end if
read(50) version_
if ((version(1).ne.version_(1)).or.(version(2).ne.version_(2)) &
 .or.(version(3).ne.version_(3))) then
  write(*,*)
  write(*,'("Warning(readdvs): different versions")')
  write(*,'(" current : ",I3.3,".",I3.3,".",I3.3)') version
  write(*,'(" file    : ",I3.3,".",I3.3,".",I3.3)') version_
  write(*,'(" in file ",A)') trim(fname)
end if
read(50) nspecies_
if (nspecies.ne.nspecies_) then
  write(*,*)
  write(*,'("Error(readdvs): differing nspecies")')
  write(*,'(" current : ",I4)') nspecies
  write(*,'(" file    : ",I4)') nspecies_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
read(50) lmmaxo_
if (lmmaxo.ne.lmmaxo_) then
  write(*,*)
  write(*,'("Error(readdvs): differing lmmaxo")')
  write(*,'(" current : ",I4)') lmmaxo
  write(*,'(" file    : ",I4)') lmmaxo_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
do js=1,nspecies
  read(50) natoms_
  if (natoms(js).ne.natoms_) then
    write(*,*)
    write(*,'("Error(readdvs): differing natoms for species ",I4)') js
    write(*,'(" current : ",I4)') natoms(js)
    write(*,'(" file    : ",I4)') natoms_
    write(*,'(" in file ",A)') trim(fname)
    write(*,*)
    stop
  end if
  read(50) nrmt_
  if (nrmt(js).ne.nrmt_) then
    write(*,*)
    write(*,'("Error(readdvs): differing nrmt for species ",I4)') js
    write(*,'(" current : ",I6)') nrmt(js)
    write(*,'(" file    : ",I6)') nrmt_
    write(*,'(" in file ",A)') trim(fname)
    write(*,*)
    stop
  end if
end do
read(50) ngridg_
if ((ngridg(1).ne.ngridg_(1)).or.(ngridg(2).ne.ngridg_(2)).or. &
 (ngridg(3).ne.ngridg_(3))) then
  write(*,*)
  write(*,'("Error(readdvs): differing ngridg")')
  write(*,'(" current : ",3I6)') ngridg
  write(*,'(" file    : ",3I6)') ngridg_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
read(50) zfmt,dvsir
do jas=1,natmtot
  js=idxis(jas)
  call zfmtpack(.true.,nrmt(js),nrmti(js),zfmt(:,:,jas),dvsmt(:,jas))
end do
close(50)
deallocate(zfmt)
return
end subroutine

