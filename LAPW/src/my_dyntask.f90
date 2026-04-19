subroutine my_dyntask(fnum, fext)
  use modmain, only: natoms, nqpt, nspecies, idxas
  use modphonon, only: ipph, isph, iaph, iqph, iasph, iq0, tphq0
  implicit none
  ! arguments
  integer, intent(in) :: fnum
  character(*), intent(out) :: fext
  ! local variables
  logical exist

  ! Some variables for phonon calculations
  ! iqph : current phonon calc q-point
  ! isph :                     species index
  ! iaph :                     atom index
  ! iasph:                     combined atom species index
  ! ipph :                     polarization index (x,y,z)

  ! search for file
  do ipph = 1,3
    do isph = 1,nspecies
      do iaph = 1,natoms(isph)
        do iqph = 1,nqpt
          ! construct the phonon file extension
          call phfext(iqph,isph,iaph,ipph,fext)
          ! determine if the DYN file with this extension exists
          write(*,*) 'In dyntask: filename = ', 'DYN'//trim(fext)
          inquire(file='DYN'//trim(fext), exist=exist)
          if (.not. exist) then
            open(fnum,file='DYN'//trim(fext),form='FORMATTED')
            iasph = idxas(iaph,isph)
            goto 10
          endif
        enddo
      enddo
    enddo
  enddo
  iqph = 0;
  isph = 0;
  iaph = 0;
  iasph = 0;
  ipph = 0
  write(*,*)
  write(*,'("Info(dyntask): nothing more to do")')
  10 continue

  if (iqph == 0) then
    fext = '.OUT'
  else
    call phfext(iqph,isph,iaph,ipph,fext)
  endif

  ! set the q=0 flag
  if(iqph == iq0) then
    tphq0 = .true.
  else
    tphq0 = .false.
  endif
  return

end subroutine

