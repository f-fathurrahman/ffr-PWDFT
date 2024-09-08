
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhouplot
use modmain
use modulr
implicit none
! initialise universal variables
call init0
call init1
! initialise the ultra long-range variables
call initulr
! read in the density from STATE_ULR.OUT
call readstulr
! write the density plot to file
select case(task)
case(731)
  open(50,file='RHOU1D.OUT',form='FORMATTED')
  open(51,file='RHOULINES.OUT',form='FORMATTED')
  call plotu1d(50,51,1,rhoqmt,rhoqir)
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(rhouplot):")')
  write(*,'(" 1D ultra long-range density plot written to RHOU1D.OUT")')
  write(*,'(" vertex location lines written to RHOULINES.OUT")')
case(732)
  open(50,file='RHOU2D.OUT',form='FORMATTED')
  call plotu2d(.false.,50,1,rhoqmt,rhoqir)
  close(50)
  write(*,*)
  write(*,'("Info(rhouplot): 2D ultra long-range density plot written to &
   &RHOU2D.OUT")')
case(733)
  open(50,file='RHOU3D.OUT',form='FORMATTED')
  call plotu3d(50,1,rhoqmt,rhoqir)
  close(50)
  write(*,*)
  write(*,'("Info(rhouplot): 3D ultra long-range density plot written to &
   &RHOU3D.OUT")')
end select
return
end subroutine

