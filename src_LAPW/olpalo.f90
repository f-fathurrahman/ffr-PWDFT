SUBROUTINE olpalo(ias,ngp,apwalm,ld,o)
  USE m_gkvectors, ONLY: ngkmax
  USE m_muffin_tins, ONLY: lmmaxapw, idxlm
  USE m_atoms, ONLY: idxis
  USE m_apwlo, ONLY: apwordmax, nlorb, apword, lorbl
  USE m_hamiltonian, ONLY: oalo, idxlo
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ias,ngp
  COMPLEX(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw)
  INTEGER, intent(in) :: ld
  COMPLEX(8), intent(inout) :: o(*)
  ! local variables
  INTEGER is,ilo,io
  INTEGER l,m,lm,i,j,k
  is=idxis(ias)
  DO ilo=1,nlorb(is)
    l=lorbl(ilo,is)
    DO m=-l,l
      lm=idxlm(l,m)
      j=ngp+idxlo(lm,ilo,ias)
      k=(j-1)*ld
      DO i=1,ngp
        k=k+1
        DO io=1,apword(l,is)
          o(k)=o(k)+conjg(apwalm(i,io,lm))*oalo(io,ilo,ias)
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 

