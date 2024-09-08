!--------------------------
subroutine initialize_elk()
!--------------------------
use modmain
use modmpi
use modomp
use modvars
implicit none
! local variables
logical exist
integer itask
! initialise MPI execution environment
call mpi_init(ierror)
! duplicate mpi_comm_world
call mpi_comm_dup(mpi_comm_world,mpicom,ierror)
! determine the number of MPI processes
call mpi_comm_size(mpicom,np_mpi,ierror)
! determine the local MPI process number
call mpi_comm_rank(mpicom,lp_mpi,ierror)
! determine if the local process is the master
if (lp_mpi.eq.0) then
  mp_mpi=.true.
  write(*,*)
  write(*,'("Elk code version ",I1.1,".",I1.1,".",I2.2," started")') version
else
  mp_mpi=.false.
end if
! read input files
call readinput
! initialise OpenMP variables
call omp_init
! initialise the MKL library
call mkl_init
! initialise the OpenBLAS library
call oblas_init
! initialise the BLIS library
call blis_init
if (mp_mpi) then
  write(*,*)
  write(*,'("Number of MPI processes : ",I6)') np_mpi
  write(*,'("Number of OpenMP threads per MPI process : ",I4)') maxthd
  write(*,'("Total number of threads : ",I6)') np_mpi*maxthd
  write(*,'("Maximum OpenMP nesting level : ",I4)') maxlvl
end if
! delete the VARIABLES.OUT file
call delvars
! write version number to VARIABLES.OUT
call writevars('version',nv=3,iva=version)
! check if Elk is already running in this directory
if (mp_mpi) then
  inquire(file='RUNNING',exist=exist)
  if (exist) then
    write(*,*)
    write(*,'("Info(elk): several copies of Elk may be running in this path")')
    write(*,'("(this could be intentional, or result from a previous crash,")')
    write(*,'(" or arise from an incorrect MPI compilation)")')
  else
    open(50,file='RUNNING')
    close(50)
  end if
end if


end subroutine