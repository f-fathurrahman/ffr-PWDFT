SUBROUTINE olpfv(nmatp,ngp,igpig,apwalm,o)
  USE m_gkvectors, ONLY: ngkmax
  USE m_muffin_tins, ONLY: lmmaxapw
  USE m_atoms, ONLY: natmtot
  USE m_apwlo, ONLY: apwordmax
  USE m_hamiltonian, ONLY: tefvr
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nmatp,ngp,igpig(ngkmax)
  COMPLEX(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  COMPLEX(8), intent(out) :: o(nmatp,nmatp)
  ! local variables
  INTEGER ias,i

  ! zero the upper triangular part of the matrix
  DO i=1,nmatp
    o(1:i,i)=0.d0
  ENDDO 
  o(:,:) = 0.0d0
  !write(*,*) 'sum o before olpaa: ', sum( abs(o(1:nmatp,1:nmatp)) )
  DO ias=1,natmtot
    CALL olpaa(tefvr,ias,ngp,apwalm(:,:,:,ias),nmatp,o)
    write(*,*) 'olpaa: ia, sum abs O = ', ias, sum( abs(o(1:nmatp,1:nmatp)) )
  ENDDO 
  write(*,*) 'Before olpistl, sum abs O: ', sum( abs(o(1:nmatp,1:nmatp)) )
  CALL olpistl(ngp,igpig,nmatp,o)  
  write(*,*) 'After olpistl, sum abs O: ', sum( abs(o(1:nmatp,1:nmatp)) )
  
  DO ias=1,natmtot
    CALL olpalo(ias,ngp,apwalm(:,:,:,ias),nmatp,o)
    CALL olplolo(ias,ngp,nmatp,o)
  ENDDO 
  
  RETURN 
END SUBROUTINE 

