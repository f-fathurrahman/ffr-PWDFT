SUBROUTINE hmlalo(ias,ngp,apwalm,ld,h)
  USE m_atoms, ONLY: idxis
  USE m_muffin_tins, ONLY: lmmaxapw, idxlm, lmaxo, lmaxapw
  USE m_apwlo, ONLY: apwordmax, nlorb, apword, lorbl
  USE m_gkvectors, ONLY: ngkmax
  USE m_hamiltonian, ONLY: hloa, gntyry, idxlo
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ias,ngp
  COMPLEX(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
  INTEGER, intent(in) :: ld
  COMPLEX(8), intent(inout) :: h(*)
  ! local variables
  INTEGER is,io,ilo
  INTEGER l1,l2,l3,m1,m2,m3
  INTEGER lm1,lm2,lm3,i,j,k
  COMPLEX(8) z1
  is=idxis(ias)
  DO ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    DO m1=-l1,l1
      lm1=idxlm(l1,m1)
      j=ngp+idxlo(lm1,ilo,ias)
      lm3=0
      DO l3=0,lmaxapw
        DO m3=-l3,l3
          lm3=lm3+1
          DO io=1,apword(l3,is)
            z1=0.d0
            DO l2=0,lmaxo
              IF(mod(l1+l2+l3,2) == 0) THEN 
                DO m2=-l2,l2
                  lm2=idxlm(l2,m2)
                  z1=z1+gntyry(lm2,lm3,lm1)*hloa(lm2,io,l3,ilo,ias)
                ENDDO 
              ENDIF 
            ENDDO 
            ! note that what is actually computed is the Hermitian conjugate of <lo|H|APW>
            IF(abs(dble(z1))+abs(aimag(z1)) > 1.d-14) THEN 
              k=(j-1)*ld
              DO i=1,ngp
                k=k+1
                h(k)=h(k)+conjg(z1*apwalm(i,io,lm3))
              ENDDO 
            ENDIF 
          ENDDO 
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO 
  RETURN
END SUBROUTINE 
