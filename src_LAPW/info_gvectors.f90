SUBROUTINE info_gvectors()
  USE m_gvectors
  USE m_gkvectors
  IMPLICIT NONE 
  INTEGER :: ig

  WRITE(*,*)
  WRITE(*,*) '---------------------------'
  WRITE(*,*) 'Info Gvectors and Gkvectors'
  WRITE(*,*) '---------------------------'

  WRITE(*,'(1x,A,F18.10)') 'rgkmax = ', rgkmax
  WRITE(*,'(1x,A,F18.10)') 'gkmax  = ', gkmax
  WRITE(*,*) 'ngk = ', ngk
  WRITE(*,'(1x,A,F18.10)') 'gmaxvr = ', gmaxvr
  WRITE(*,*)
  WRITE(*,*) 'ngvec = ', ngvec
  WRITE(*,*) 'ngtot = ', ngtot 
  WRITE(*,*) 'size(gc) = ', size(gc)
  WRITE(*,'(1x,A,3I8)') 'ngridg = ', ngridg(:)
  WRITE(*,'(1x,A,3I8)') 'ngdc   = ', ngdc(:)

  WRITE(*,*)
  WRITE(*,*) 'Some gvectors less than gmaxvr: '
  DO ig = 1,5
    WRITE(*,'(1x,I8,3F18.10,A,2F18.10)') ig, vgc(:,ig), '|', gc(ig), gc(ig)**2
  ENDDO 
  WRITE(*,*) '....'
  DO ig = ngvec-5,ngvec
    WRITE(*,'(1x,I8,3F18.10,A,2F18.10)') ig, vgc(:,ig), '|', gc(ig), gc(ig)**2
  ENDDO 
  WRITE(*,*)

END SUBROUTINE 
