
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modomp

! maximum number of OpenMP threads available
integer maxthd
! maximum number of OpenMP threads for the first nesting level
integer maxthd1
! maximum number of threads available to MKL
integer maxthdmkl
! maximum OpenMP nesting level
integer maxlvl
! number of active OpenMP threads for each nesting level
integer, allocatable :: nathd(:)

interface

integer function omp_get_num_procs()
end function

integer function omp_get_max_threads()
end function

integer function omp_get_num_threads()
end function

integer function omp_get_thread_num()
end function

logical function omp_get_nested()
end function

integer function omp_get_max_active_levels()
end function

logical function omp_get_dynamic()
end function

integer function omp_get_level()
end function

subroutine omp_set_num_threads(num_threads)
integer, intent(in) :: num_threads
end subroutine

subroutine omp_set_nested(nested)
logical, intent(in) :: nested
end subroutine

subroutine omp_set_max_active_levels(max_levels)
integer, intent(in) :: max_levels
end subroutine

subroutine omp_set_dynamic(dynamic_threads)
logical, intent(in) :: dynamic_threads
end subroutine

end interface

contains

subroutine omp_init
implicit none
if (maxthd.lt.0) then
! set the number of threads equal to the number of processors
  maxthd=omp_get_num_procs()
  call omp_set_num_threads(maxthd)
else if (maxthd.eq.0) then
! use the system default number of threads
  maxthd=omp_get_max_threads()
else
! use the number of threads specified in the input file
  call omp_set_num_threads(maxthd)
end if
if (maxthd1.le.0) then
  maxthd1=maxthd
else
  maxthd1=min(maxthd1,maxthd)
end if
! switch off dynamic allocation of threads
call omp_set_dynamic(.false.)
! allow nested parallelism
call omp_set_nested(.true.)
! set the maximum nesting level
call omp_set_max_active_levels(maxlvl)
! allocate the number of active threads array
if (allocated(nathd)) deallocate(nathd)
allocate(nathd(0:maxlvl))
! initialise the number of active threads
call omp_reset
return
end subroutine

subroutine omp_reset
implicit none
! number of active threads at each nesting level
nathd(0)=1
nathd(1:)=0
return
end subroutine

subroutine holdthd(nloop,nthd)
implicit none
! arguments
integer, intent(in) :: nloop
integer, intent(out) :: nthd
! local variables
integer lvl,na,n
! current nesting level
lvl=omp_get_level()
if ((lvl.lt.0).or.(lvl.ge.maxlvl)) then
  nthd=1
  return
end if
!$OMP CRITICAL(holdthd_)
! determine number of active threads at the current nesting level
na=nathd(lvl)
na=max(min(na,maxthd),1)
! number of threads allowed for this loop
nthd=maxthd/na
if (mod(maxthd,na).gt.0) nthd=nthd+1
if (lvl.eq.0) nthd=min(nthd,maxthd1)
nthd=max(min(nthd,maxthd,nloop),1)
! add to number of active threads in next nesting level
n=nathd(lvl+1)+nthd
n=max(min(n,maxthd),0)
nathd(lvl+1)=n
!$OMP END CRITICAL(holdthd_)
return
end subroutine

subroutine freethd(nthd)
implicit none
! arguments
integer, intent(in) :: nthd
! local variables
integer lvl,n
! current nesting level
lvl=omp_get_level()
if ((lvl.lt.0).or.(lvl.ge.maxlvl)) return
!$OMP CRITICAL(freethd_)
! subtract from the number of active threads in next nesting level
n=nathd(lvl+1)-nthd
n=max(min(n,maxthd),0)
nathd(lvl+1)=n
!$OMP END CRITICAL(freethd_)
return
end subroutine

end module

