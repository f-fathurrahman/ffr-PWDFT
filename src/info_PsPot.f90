SUBROUTINE info_PsPot()
  USE m_atoms, ONLY : Nspecies, SpeciesSymbols
  USE m_PsPot, ONLY : Ps => Ps_HGH_Params, &
                      PsPot_FilePath
  use m_Ps_HGH, ONLY : info_Ps_HGH_Params
  IMPLICIT NONE 
  INTEGER :: isp

  WRITE(*,*)
  WRITE(*,'(9x,A)') 'Pseudopotential Information'
  WRITE(*,'(9x,A)') '---------------------------'
  DO isp = 1,Nspecies
    WRITE(*,*)
    WRITE(*,*) 'Species: ', trim(SpeciesSymbols(isp))
    WRITE(*,*) 'Pseudopotential file: ', trim(PsPot_FilePath(isp))
    WRITE(*,*)
    CALL info_Ps_HGH_Params( Ps(isp) )
  ENDDO 

END SUBROUTINE 
