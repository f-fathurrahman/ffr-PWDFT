
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mkl_init
use modomp
implicit none
! set the initial number of MKL threads equal to one
call mkl_set_num_threads(1)
! set the maximum number of threads available to MKL
if (maxthdmkl.le.0) maxthdmkl=maxthd
! enable dynamic thread allocation
call mkl_set_dynamic(.true.)
return
end subroutine

