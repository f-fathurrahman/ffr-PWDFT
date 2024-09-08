
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writephn
use modmain
use modphonon
implicit none
! local variables
integer iq,i,j,is,ia,ip
! allocatable arrays
real(8), allocatable :: w(:)
complex(8), allocatable :: dynq(:,:,:),dynr(:,:,:)
complex(8), allocatable :: dynp(:,:),ev(:,:)
! initialise universal variables
call init0
call init2
allocate(w(nbph))
allocate(dynq(nbph,nbph,nqpt))
allocate(dynr(nbph,nbph,nqptnr))
allocate(dynp(nbph,nbph))
allocate(ev(nbph,nbph))
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! Fourier transform the dynamical matrices to real-space
call dynqtor(dynq,dynr)
open(50,file='PHONON.OUT',form='FORMATTED')
do iq=1,nphwrt
  call dynrtoq(vqlwrt(:,iq),dynr,dynp)
  call dynev(dynp,w,ev)
  write(50,*)
  write(50,'(I6,3G18.10," : q-point, vqlwrt")') iq,vqlwrt(:,iq)
  do j=1,nbph
    write(50,*)
    write(50,'(I6,G18.10," : mode, frequency")') j,w(j)
    i=0
    do is=1,nspecies
      do ia=1,natoms(is)
        do ip=1,3
          i=i+1
          if (i.eq.1) then
            write(50,'(3I4,2G18.10," : species, atom, polarisation, &
             &eigenvector")') is,ia,ip,ev(i,j)
          else
            write(50,'(3I4,2G18.10)') is,ia,ip,ev(i,j)
          end if
        end do
      end do
    end do
  end do
  write(50,*)
end do
close(50)
write(*,*)
write(*,'("Info(writephn): phonon frequencies and eigenvectors written to &
 &PHONON.OUT")')
write(*,'(" for all q-vectors in the phwrite list")')
deallocate(w,dynq,dynr,dynp,ev)
return
end subroutine

