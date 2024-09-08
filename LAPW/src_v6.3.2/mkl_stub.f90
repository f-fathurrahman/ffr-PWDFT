
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Stub routines for Intel MKL

subroutine mkl_set_num_threads(num_threads)
implicit none
integer, intent(in) :: num_threads
return
end subroutine

subroutine mkl_set_dynamic(dynamic)
implicit none
logical, intent(in) :: dynamic
return
end subroutine

