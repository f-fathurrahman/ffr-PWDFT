SUBROUTINE add_lorb_cnd()
  USE m_atoms, ONLY: nspecies
  USE m_apwlo, ONLY: lorbcnd, lorbve, lorbdm, nlorb, lorbe0, lorbord, lorbl, &
               maxlorb, lorbordc
  USE m_muffin_tins, ONLY: lmaxo
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is, nlo, l, io

  IF(.not. lorbcnd) RETURN

  ! add conduction local-orbitals to each species
  DO is = 1,nspecies
    nlo = nlorb(is)
    DO l = 0,lmaxo
      nlo = nlo + 1
      IF(nlo > maxlorb) THEN 
        WRITE(*,*)
        WRITE(*,'("Error(addlorbcnd): nlorb too large : ",I8)') nlo
        WRITE(*,'(" for species ",I4)') is
        WRITE(*,'("Adjust maxlorb in modmain and recompile code")')
        WRITE(*,*)
        STOP 
      ENDIF 
      lorbl(nlo,is)=l
      lorbord(nlo,is)=lorbordc
      DO io = 1,lorbordc
        lorbe0(io,nlo,is) = 0.15d0
        lorbdm(io,nlo,is) = io - 1
        lorbve(io,nlo,is) = .true.
      ENDDO 
    ENDDO 
    nlorb(is) = nlo
  ENDDO 
  RETURN 
END SUBROUTINE 