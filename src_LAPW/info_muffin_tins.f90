subroutine info_muffin_tins()
  USE m_atoms
  USE m_muffin_tins

  write(*,*)
  write(*,*) '----------------'  
  write(*,*) 'Muffin Tins Info'
  write(*,*) '----------------'
  write(*,*)

  write(*,*) 'Muffin Tin radius'
  do isp=1,Nspecies
    write(*,'(1x,I4,F18.10)') isp, rmt(isp)
  enddo

  write(*,*)
  write(*,*) 'Angular-momentum'
  write(*,'(1x,A,I4)') 'lmaxapw  = ', lmaxapw
  write(*,'(1x,A,I4)') 'lmmaxapw = ', lmmaxapw
  write(*,'(1x,A,I4)') 'lmaxo    = ', lmaxo
  write(*,'(1x,A,I4)') 'lmmaxo   = ', lmmaxo
  write(*,'(1x,A,I4)') 'lmaxi    = ', lmaxi
  write(*,'(1x,A,I4)') 'lmmaxi   = ', lmmaxi
end subroutine