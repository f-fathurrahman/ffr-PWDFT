subroutine info_muffin_tins()
  USE m_atoms
  USE m_muffin_tins

  write(*,*)
  write(*,*) '---------------------'  
  write(*,*) 'Muffin Tins (MT) Info'
  write(*,*) '---------------------'
  write(*,*)

  write(*,*) 'rmt (MT radius)'
  do isp=1,Nspecies
    write(*,'(1x,I4,F18.10)') isp, rmt(isp)
  enddo

  write(*,*)
  write(*,*) 'nrmt, nrmti (fine MT radial points, outer inner)'
  do isp=1,Nspecies
    write(*,'(1x,I4,2I8)') isp, nrmt(isp), nrmti(isp)
  enddo

  write(*,*)
  write(*,*) 'nrcmt, nrcmti (coarse MT radial points)'
  do isp=1,Nspecies
    write(*,'(1x,I4,2I8)') isp, nrcmt(isp), nrcmti(isp)
  enddo

  write(*,*)
  write(*,*) 'npmt, npmti (packed fine MT radial points)'
  do isp=1,Nspecies
    write(*,'(1x,I4,2I8)') isp, npmt(isp), npmti(isp)
  enddo

  write(*,*)
  write(*,*) 'npcmt, npmcti (packed coarse MT radial points)'
  do isp=1,Nspecies
    write(*,'(1x,I4,2I8)') isp, npcmt(isp), npcmti(isp)
  enddo


  write(*,*)
  write(*,*) 'lradstp  = ', lradstp
  write(*,*) 'nrmtmax  = ', nrmtmax
  write(*,*) 'nrcmtmax = ', nrcmtmax
  write(*,'(1x,A,F18.10)') 'omegamt  = ', omegamt
  write(*,'(1x,A,F18.10)') 'rmtdelta = ', rmtdelta

  write(*,*)
  write(*,*) 'Angular-momentum'
  write(*,'(1x,A,I4)') 'lmaxapw  = ', lmaxapw
  write(*,'(1x,A,I4)') 'lmmaxapw = ', lmmaxapw
  write(*,'(1x,A,I4)') 'lmaxo    = ', lmaxo
  write(*,'(1x,A,I4)') 'lmmaxo   = ', lmmaxo
  write(*,'(1x,A,I4)') 'lmaxi    = ', lmaxi
  write(*,'(1x,A,I4)') 'lmmaxi   = ', lmmaxi


  write(*,*) 'Some rlmt:'
  write(*,*) 'rlmt = '
  write(*,*) 2, 0, 1, rlmt(2,0,1)
  write(*,*) 2, -1, 1, rlmt(2,-1,1)
  write(*,*) 2, -7, 1, rlmt(2,-7,1)
end subroutine