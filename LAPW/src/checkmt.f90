!   Checks for muffin-tins which are too close together or intersecting. If any
!   such muffin-tins are found THEN  the radii of their associated atomic species
!   are adjusted so that the minimum distance between their surfaces is
!   {\tt rmtdelta}.
SUBROUTINE checkmt()
  USE m_lattice, ONLY: epslat
  USE m_atoms, ONLY: nspecies
  USE m_atomic_species, ONLY: spsymb
  USE m_muffin_tins, ONLY: rmt, rmtdelta
  !
  IMPLICIT NONE 
  ! local variables
  INTEGER is,js
  REAL(8) dmin,t1,t2
  ! automatic arrays
  REAL(8) rmt0(nspecies)
  
  rmt0(1:nspecies)=rmt(1:nspecies)
  
  10 CONTINUE 
  
  ! find the minimum distance between muffin-tin surfaces
  CALL mtdmin(is, js, dmin)
  
  ! adjust muffin-tin radii if required
  IF(dmin < rmtdelta-epslat) THEN 
    WRITE(*,*) 'INFO: Adjusting rmt'
    t1 = rmt(is) + rmt(js)
    t2 = (t1 + dmin - rmtdelta)/t1
    rmt(is) = rmt(is)*t2
    IF(is /= js) rmt(js) = rmt(js)*t2
    GOTO 10
  ENDIF 
  
  DO is=1,nspecies
    IF(rmt(is) < 0.25d0) THEN 
      WRITE(*,*)
      WRITE(*,'(1x,"ERROR: muffin-tin radius too small for species ",I4," (",A,")")') &
        is, trim(spsymb(is))
      WRITE(*,'(" Radius : ",G18.10)') rmt(is)
      WRITE(*,*)
      STOP 
    ENDIF 
    IF(rmt(is) < rmt0(is)) THEN 
      WRITE(*,'(1x,"INFO: reduced muffin-tin radius of species ",I3," (",A,") from ",F8.4," to ",F8.4)') &
        is, trim(spsymb(is)), rmt0(is), rmt(is)
    ENDIF 
  ENDDO 
  RETURN 
END SUBROUTINE 
