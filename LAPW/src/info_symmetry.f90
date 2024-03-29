SUBROUTINE info_symmetry()
  USE m_symmetry
  USE m_atoms, ONLY: natmtot
  IMPLICIT NONE
  INTEGER :: ia

  WRITE(*,*)
  WRITE(*,*) '---------------------------'
  WRITE(*,*) 'Info symmetry'
  WRITE(*,*) '---------------------------'

  WRITE(*,*)
  WRITE(*,'(1x,A,I4)') 'nsymlat  = ', nsymlat
  WRITE(*,'(1x,A,I4)') 'nsymcrys = ', nsymcrys
  WRITE(*,*)
  DO ia=1,natmtot
    WRITE(*,'(1x,A,I4,A,I4)') 'Site: ', ia, ' nsymsite: ', nsymsite(ia)
  ENDDO

END SUBROUTINE 
