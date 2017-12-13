SUBROUTINE info_PWGrid()

  USE m_PWGrid
  USE m_cell
  USE m_realspace
  IMPLICIT NONE 
  INTEGER :: i
  
  WRITE(*,*)
  WRITE(*,'(1x,A)') 'Unit cell vectors:'
  DO i = 1,3
    WRITE(*,'(1x,3F18.10)') LatVecs(i,:)
  ENDDO 

  WRITE(*,*)
  WRITE(*,'(1x,A)') 'Unit reciprocal cell vectors:'
  DO i = 1,3
    WRITE(*,'(1x,3F18.10)') RecVecs(i,:)
  ENDDO 

  WRITE(*,*)
  WRITE(*,'(1x,A,F18.10)') 'ecutwfc (Ha)     = ', ecutwfc
  WRITE(*,'(1x,A,F18.10)') 'ecutrho (Ha)     = ', ecutrho
  WRITE(*,'(1x,A,F18.10)') 'Unit cell volume = ', CellVolume
  WRITE(*,'(1x,A,3I8)')    'Sampling point   = ', Ns(:)
  WRITE(*,'(1x,A,F18.10)') 'dVol             = ', dVol
  WRITE(*,'(1x,A,I8)')     'Npoints          = ', Npoints
  WRITE(*,'(1x,A,I8)')     'Ngwx             = ', Ngwx

END SUBROUTINE 

