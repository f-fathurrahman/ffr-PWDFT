SUBROUTINE info_crystal()
  USE m_lattice
  USE m_constants, ONLY: pi
  USE m_atoms
  IMPLICIT NONE 
  INTEGER :: i, is, ia
  REAL(8) :: m(3,3)

  WRITE(*,*)
  WRITE(*,*) '---------------------------'
  WRITE(*,*) 'Info crystal'
  WRITE(*,*) '---------------------------'

  WRITE(*,*)
  WRITE(*,*) 'avec matrix'
  DO i=1,3
    WRITE(*,'(1x,3F18.10)') avec(i,:)
  ENDDO 

  WRITE(*,*)
  WRITE(*,*) 'bvec matrix'
  DO i=1,3
    WRITE(*,'(1x,3F18.10)') bvec(i,:)
  ENDDO 

  m = matmul(avec, transpose(bvec))/(2.d0*pi)
  WRITE(*,*)
  WRITE(*,*) 'avec*bvec matrix = '
  DO i=1,3
    WRITE(*,'(1x,3F18.10)') m(i,:)
  ENDDO 

  WRITE(*,*)
  WRITE(*,*) 'atposl: '
  DO is = 1,nspecies
    DO ia = 1,natoms(is)
      WRITE(*,'(1x,3F18.10)') atposl(:,ia,is)
    ENDDO 
  ENDDO 

  WRITE(*,*)
  WRITE(*,*) 'atposc: '
  DO is = 1,nspecies
    DO ia = 1,natoms(is)
      WRITE(*,'(1x,3F18.10)') atposc(:,ia,is)
    ENDDO 
  ENDDO 

END SUBROUTINE 
