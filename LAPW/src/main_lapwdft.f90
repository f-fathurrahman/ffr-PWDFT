PROGRAM lapwdft
  use modmain, only: tasks, ntasks, task
  integer :: itask

  CALL read_input()

  write(*,*) 'tasks = ', tasks

  ! perform the tasks
  do itask = 1,ntasks
    task = tasks(itask) ! this will set global variable task
    !
    select case(task)
    case(0,1)
      call gndstate()
    case(120)
      call writepmat()
    case(121)
      call dielectric()
    case(205)
      call phonon()
    case(210)
      call phdos()
    case(220)
      call phdisp()
    case default
      write(*,*)
      write(*,'("Error(elk): task not defined or not supported: ",I8)') task
      write(*,*)
      stop
    end select
  enddo

END PROGRAM 

