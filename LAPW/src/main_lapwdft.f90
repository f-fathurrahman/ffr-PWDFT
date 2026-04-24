PROGRAM lapwdft
  use modmain, only: tasks, ntasks, task
  integer :: itask

  CALL read_input()

  write(*,*) 'tasks = ', tasks(1:ntasks)

  ! perform the tasks
  do itask = 1,ntasks
    task = tasks(itask) ! this will set global variable task
    !
    select case(task)
    case(0,1)
      call gndstate()
    case(2,3)
      call geomopt()
    case(5)
      call hartfock()
    case(10)
      call writedos()
    case(20,21,22,23)
      call bandstr()
    case(25)
      call effmass()
    case(120)
      call my_writepmat()
    case(121)
      call my_dielectric()
    case(180)
      call writeepsinv()
    case(185)
      call writehmlbse()
    case(186)
      call writeevbse()
    case(187)
      call dielectric_bse()
    case(205)
      !call phonon()
      call my_phonon()
    case(210)
      call phdos()
    case(220)
      call phdisp()
    !case(300)
    !  call rdmft()
    case(320)
      call tddftlr()
    !case(330,331)
    !  call tddftsplr()
    !case(341,342,343)
    !  call wxcplot()
    case default
      write(*,*)
      write(*,'("Error(elk): task not defined or not supported: ",I8)') task
      write(*,*)
      stop
    end select
  enddo

END PROGRAM 

