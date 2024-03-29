SUBROUTINE genffacgp(is,gpc,ffacgp)
  USE m_constants, ONLY: fourpi
  USE m_gvectors, ONLY: ngtot
  USE m_lattice, ONLY: epslat, omega
  USE m_muffin_tins,  ONLY: rmt
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: is
  REAL(8), intent(in) :: gpc(ngtot)
  REAL(8), intent(out) :: ffacgp(ngtot)
  ! local variables
  INTEGER ig
  REAL(8) t1,t2
  t1=fourpi/omega
  DO ig=1,ngtot
    IF(gpc(ig).gt.epslat) THEN 
      t2=gpc(ig)*rmt(is)
      ffacgp(ig)=t1*(sin(t2)-t2*cos(t2))/(gpc(ig)**3)
    else
      ffacgp(ig)=(t1/3.d0)*rmt(is)**3
    ENDIF 
  ENDDO 
  RETURN 
END SUBROUTINE 

