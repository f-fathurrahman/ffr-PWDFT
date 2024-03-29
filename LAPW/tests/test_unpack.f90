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

  CALL test_unpack()

END PROGRAM

!-----------------------
subroutine test_unpack()
!-----------------------
  USE m_density_pot_xc, ONLY: rhomt
  USE m_atoms, ONLY: natmtot, idxis
  USE m_muffin_tins, ONLY: nrmtmax, lmmaxo, nrmt, nrmti
  implicit none
  integer :: is, ias, ir, lm
  real(8), allocatable :: rfmt1(:,:,:)

  ! unpack the muffin-tin function
  ALLOCATE(rfmt1(lmmaxo,nrmtmax,natmtot))
  DO ias = 1,natmtot
    is = idxis(ias)
    CALL my_rfmtpack(.false., nrmt(is), nrmti(is), rhomt(:,ias), rfmt1(:,:,ias))
  ENDDO 

  write(*,*) 'lmmaxo = '
  write(*,*) 'shape rhomt = ', shape(rhomt)
  write(*,*) 'shape rfmt1 = ', shape(rfmt1)

  do lm = 1,lmmaxo
    write(1000,*) 'lm = ', lm
    ! Fix ias=1
    do ir = 1,nrmt(1)
      write(1000,*) ir, rfmt1(lm,ir,1)
    enddo
  enddo

  write(*,*) 'Pass here ...'

  deallocate(rfmt1)

end subroutine

include 'my_rfmtpack.f90'