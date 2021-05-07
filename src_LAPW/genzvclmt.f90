SUBROUTINE genzvclmt(nr,nri,ld1,rl,wpr,ld2,zrhomt,zvclmt)
  USE m_atoms, ONLY: nspecies, natmtot, idxis
  USE m_muffin_tins, ONLY: lmaxo
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr(nspecies),nri(nspecies)
  INTEGER, intent(in) :: ld1
  REAL(8), intent(in) :: rl(ld1,-lmaxo-1:lmaxo+2,nspecies)
  REAL(8), intent(in) :: wpr(4,ld1,nspecies)
  INTEGER, intent(in) :: ld2
  COMPLEX(8), intent(in) :: zrhomt(ld2,natmtot)
  COMPLEX(8), intent(out) :: zvclmt(ld2,natmtot)
  ! local variables
  INTEGER :: is,ias
  DO ias=1,natmtot
    is=idxis(ias)
    CALL zpotclmt(nr(is),nri(is),ld1,rl(:,:,is),wpr(:,:,is),zrhomt(:,ias), &
     zvclmt(:,ias))
  ENDDO 
  RETURN 
END SUBROUTINE 
