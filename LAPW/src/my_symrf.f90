SUBROUTINE my_symrf(nr, nri, np, ld, rfmt, rfir)
  ! !USES:
  USE m_atoms, ONLY: nspecies, natmtot
  USE m_gvectors, ONLY: ngtot
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
  INTEGER, intent(in) :: ld
  REAL(8), intent(inout) :: rfmt(ld,natmtot), rfir(ngtot)
  
  write(*,*)
  write(*,*) '<div> ENTER my_symrf'
  write(*,*)

  write(*,*) 'ld (should be equal to npmtmax) = ', ld
  ! local variables
  CALL my_symrfmt(nr, nri, np, ld, rfmt)
  CALL my_symrfir(rfir)

  write(*,*)
  write(*,*) '</div> EXIT my_symrf'
  write(*,*)

  RETURN 
END SUBROUTINE 
