SUBROUTINE my_symrvf(tspin,tnc,nr,nri,np,ld,rvfmt,rvfir)
  ! !USES:
  USE modmain
  IMPLICIT NONE 
  ! arguments
  logical, intent(in) :: tspin,tnc
  INTEGER, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
  INTEGER, intent(in) :: ld
  REAL(8), intent(inout) :: rvfmt(ld,natmtot,*),rvfir(ngtot,*)
  
  ! local variables
  CALL my_symrvfmt(tspin,tnc,nr,nri,np,ld,rvfmt)
  !write(*,*)
  !write(*,*) '>>>> EARLY RETURN in my_symrvf after my_symrvfmt'
  !RETURN
  CALL my_symrvfir(tspin,tnc,rvfir)
  RETURN 
END SUBROUTINE 
