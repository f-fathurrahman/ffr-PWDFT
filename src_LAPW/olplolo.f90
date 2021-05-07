SUBROUTINE olplolo(ias,ngp,ld,o)
  USE m_atoms, ONLY: idxis
  USE m_apwlo, ONLY: nlorb, lorbl
  USE m_muffin_tins, ONLY: idxlm
  USE m_hamiltonian, ONLY: ololo, idxlo
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ias,ngp,ld
  COMPLEX(8), intent(inout) :: o(ld,*)
  ! local variables
  INTEGER is,ilo,jlo
  INTEGER l,m,lm,i,j
  is=idxis(ias)
  DO ilo=1,nlorb(is)
    l=lorbl(ilo,is)
    DO jlo=1,nlorb(is)
      IF(lorbl(jlo,is) == l) THEN 
        DO m=-l,l
          lm=idxlm(l,m)
          i=ngp+idxlo(lm,ilo,ias)
          j=ngp+idxlo(lm,jlo,ias)
          IF(i.le.j) THEN 
            o(i,j)=o(i,j)+ololo(ilo,jlo,ias)
          ENDIF 
        ENDDO 
      ENDIF 
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 

