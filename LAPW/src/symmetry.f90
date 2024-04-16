SUBROUTINE symmetry()

  USE m_lattice, ONLY: avec, ainv
  USE m_symmetry, ONLY: tsyminv, symtype, nsymlat
  USE m_hamiltonian, ONLY: tefvr

  IMPLICIT NONE 

  WRITE(*,*) 'Setting up symmetry'

  ! inverse of the lattice vector matrix
  CALL r3minv(avec,ainv) ! ffr: CALL this again?

  ! find Bravais lattice symmetries
  CALL findsymlat()

  ! use only the identity if required
  IF(symtype == 0) nsymlat = 1

  ! find the crystal symmetries and shift atomic positions if required
  CALL my_findsymcrys()

  ! find the site symmetries
  CALL findsymsite()

  ! check if real symmetric first-variational eigen solver can be used
  IF(.not. tsyminv) tefvr = .false.

  RETURN 
END SUBROUTINE 

