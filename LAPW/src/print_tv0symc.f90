SUBROUTINE print_tv0symc()
  USE m_symmetry, only: tv0symc, nsymcrys
  IMPLICIT NONE 
  INTEGER :: i
  !
  DO i = 1,nsymcrys
    write(*,'(1x,I4,1x,L)') i, tv0symc(i)
  ENDDO 
END SUBROUTINE 

