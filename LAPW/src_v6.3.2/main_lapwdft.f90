PROGRAM lapwdft
  use modmain, only: tasks, ntasks
  integer :: itask, task

  CALL readinput()

  !CALL init0()
  !CALL init1()
  !CALL info_lattice()
  !CALL info_apwlo()

  ! perform the tasks
  do itask = 1,ntasks
    task = tasks(itask)
    !
    select case(task)
    case(0,1)
      call gndstate()
    case(120)
      call writepmat()
    case(121)
      call dielectric()
    case default
      write(*,*)
      write(*,'("Error(elk): task not defined or not supported: ",I8)') task
      write(*,*)
      stop
    end select
  enddo

END PROGRAM 

