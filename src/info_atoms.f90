SUBROUTINE info_atoms()

  USE m_atoms
  IMPLICIT NONE 
  INTEGER :: ia, isp


  WRITE(*,*)
  WRITE(*,'(9x,A)') 'Atoms Information'
  WRITE(*,'(9x,A)') '-----------------'
  WRITE(*,*)
  WRITE(*,*) 'Natoms   = ', Natoms
  WRITE(*,*) 'Nspecies = ', Nspecies
  WRITE(*,*)
  WRITE(*,*) 'Valence information:'
  DO isp = 1, Nspecies
    WRITE(*,'(1x,A5,F10.3)') adjustl(SpeciesSymbols(isp)), AtomicValences(isp)
  ENDDO

  WRITE(*,*)
  WRITE(*,*) 'Atomic coordinates: (in bohr)'
  WRITE(*,*)
  DO ia = 1,Natoms
    isp = atm2species(ia)
    WRITE(*,'(1x,A5,3F18.10)') adjustl(SpeciesSymbols(isp)), AtomicCoords(1:3,ia)
  ENDDO 

END SUBROUTINE 
