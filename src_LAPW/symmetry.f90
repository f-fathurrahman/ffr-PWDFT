subroutine symmetry

use modmain, only: avec, lsplsymc, lspnsymc, ainv
use modmain, only: vtlsymc, tsyminv, tefvr, symtype, symlat, nsymlat, nsymcrys

implicit none

! inverse of the lattice vector matrix
call r3minv(avec,ainv) ! ffr: call this again?

! find Bravais lattice symmetries
call findsymlat()

! use only the identity if required
if (symtype.eq.0) nsymlat=1

! find the crystal symmetries and shift atomic positions if required
call findsymcrys()

! find the site symmetries
call findsymsite()

  ! check if fixed spin moments are invariant under the symmetry group
  !call checkfsm()

! check if real symmetric first-variational eigen solver can be used
if (.not.tsyminv) tefvr=.false.

return
end subroutine

