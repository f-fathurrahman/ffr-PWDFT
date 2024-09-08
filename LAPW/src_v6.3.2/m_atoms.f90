MODULE m_atoms

!--------------------------!
!     atomic variables     !
!--------------------------!
! maximum allowed species
INTEGER, parameter :: maxspecies=8
! maximum allowed atoms per species
INTEGER, parameter :: maxatoms=200
! number of species
INTEGER nspecies
! number of atoms for each species
INTEGER natoms(maxspecies)
! maximum number of atoms over all the species
INTEGER natmmax
! total number of atoms
INTEGER natmtot
! index to atoms and species
INTEGER idxas(maxatoms,maxspecies)
! inverse atoms and species indices
INTEGER idxis(maxatoms*maxspecies)
INTEGER idxia(maxatoms*maxspecies)
! molecule is .true. is the system is an isolated molecule
LOGICAL molecule
! primcell is .true. if primitive unit cell is to be found automatically
LOGICAL primcell
! atomic positions in lattice coordinates
REAL(8) atposl(3,maxatoms,maxspecies)
! atomic positions in Cartesian coordinates
REAL(8) atposc(3,maxatoms,maxspecies)
! magnitude of random displacements added to the atomic positions
REAL(8) rndatposc

END MODULE 

