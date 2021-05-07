SUBROUTINE hmlfv(nmatp,ngp,igpig,vgpc,apwalm,h)
  USE m_gkvectors, ONLY: ngkmax
  USE m_muffin_tins, ONLY: lmmaxapw
  USE m_atoms, ONLY: natmtot
  USE m_apwlo, ONLY: apwordmax
  USE m_hamiltonian, ONLY: tefvr
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nmatp,ngp,igpig(ngkmax)
  REAL(8), intent(in) :: vgpc(3,ngkmax)
  COMPLEX(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  COMPLEX(8), intent(out) :: h(nmatp,nmatp)
  ! local variables
  INTEGER ias,i
  
  ! zero the upper triangular part of the matrix
  DO i=1,nmatp
    h(1:i,i) = 0.d0
  ENDDO 
  
  DO ias=1,natmtot
    CALL hmlaa(tefvr,ias,ngp,apwalm(:,:,:,ias),nmatp,h)
  ENDDO 
  CALL hmlistl(ngp,igpig,vgpc,nmatp,h)
  
  DO ias=1,natmtot
    CALL hmlalo(ias,ngp,apwalm(:,:,:,ias),nmatp,h)
    CALL hmllolo(ias,ngp,nmatp,h)
  ENDDO 
  
  RETURN 
END SUBROUTINE 
