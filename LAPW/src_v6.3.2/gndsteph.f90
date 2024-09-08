
! Copyright (C) 2019 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine gndsteph
use modmain
use modmpi
use modomp
implicit none
! local variables
integer ik

! initialise universal variables
call init0
call init1
call readstate
call genvsig
call gencore
call readfermi
call linengy
call genapwlofr
call gensocfr
call genevfsv
call occupy
! begin the self-consistent loop
if (mp_mpi) then
! open EPH_INFO.OUT file
  open(60,file='EPH_INFO.OUT',form='FORMATTED')
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop started |")')
  write(60,'("+------------------------------+")')
end if
do iscl=1,maxscl
  if (mp_mpi) then
    write(60,*)
    write(60,'("+--------------------+")')
    write(60,'("| Loop number : ",I4," |")') iscl
    write(60,'("+--------------------+")')
    flush(60)
    write(*,*)
    write(*,'("Info(gndsteph): self-consistent loop number : ",I4)') iscl
  end if

! loop over k-points
  do ik=1,nkpt
    call eveqneph(ik)

  end do


  if (iscl.ge.maxscl) then
    if (mp_mpi) then
      write(60,*)
      write(60,'("Reached self-consistent loops maximum")')
    end if
    tlast=.true.
  end if
! reset the OpenMP thread variables
  call omp_reset
end do
10 continue
if (mp_mpi) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop stopped |")')
  write(60,'("+------------------------------+")')
! close the EPH_INFO.OUT file
  close(60)
end if
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
return
end subroutine

