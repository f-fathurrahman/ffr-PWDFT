SUBROUTINE hmllolo(ias,ngp,ld,h)
  USE m_atoms, ONLY: idxis
  USE m_muffin_tins, ONLY: idxlm, lmaxo
  USE m_apwlo, ONLY: nlorb, lorbl
  USE m_hamiltonian, ONLY: gntyry, idxlo, hlolo
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ias,ngp,ld
  COMPLEX(8), intent(inout) :: h(ld,*)
  ! local variables
  INTEGER :: is,ilo,jlo
  INTEGER :: l1,l2,l3,m1,m2,m3
  INTEGER :: lm1,lm2,lm3,i,j
  COMPLEX(8) :: z1
  is=idxis(ias)
  DO jlo=1,nlorb(is)
    l3=lorbl(jlo,is)
    DO m3=-l3,l3
      lm3=idxlm(l3,m3)
      j=ngp+idxlo(lm3,jlo,ias)
      DO ilo=1,nlorb(is)
        l1=lorbl(ilo,is)
        DO m1=-l1,l1
          lm1=idxlm(l1,m1)
          i = ngp + idxlo(lm1,ilo,ias)
          IF(i <= j) THEN 
            z1=0.d0
            DO l2=0,lmaxo
              IF(mod(l1+l2+l3,2) == 0) THEN 
                DO m2=-l2,l2
                  lm2=idxlm(l2,m2)
                  z1=z1+gntyry(lm2,lm3,lm1)*hlolo(lm2,jlo,ilo,ias)
                ENDDO 
              ENDIF 
            ENDDO 
            h(i,j)=h(i,j)+z1
          ENDIF 
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 

