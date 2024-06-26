PROGRAM main
  IMPLICIT NONE 

  ! read input and initialize some variables
  CALL read_input()
  CALL init0()
  CALL init1()
  
  CALL info_crystal()
  CALL info_symmetry()
  CALL writesym()
  CALL info_gvectors()
  CALL info_muffin_tins()

  CALL rhoinit()

  CALL test_plot3d()

END PROGRAM

!-----------------------
subroutine test_plot3d()
!-----------------------
  USE m_density_pot_xc, ONLY: rhoir, rhomt
  implicit none

  call my_plot3d(1001, 1, rhomt, rhoir)
  write(*,*) 'File saved to unit 1001'

end subroutine

include 'my_plot3d.f90'
include 'my_plotpt3d.f90'
include 'my_rfplot.f90'
include 'my_rfip.f90'