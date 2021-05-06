module m_atoms

!--------------------------!
!     atomic variables     !
!--------------------------!
! maximum allowed species
integer, parameter :: maxspecies=8
! maximum allowed atoms per species
integer, parameter :: maxatoms=200
! number of species
integer nspecies
! number of atoms for each species
integer natoms(maxspecies)
! maximum number of atoms over all the species
integer natmmax
! total number of atoms
integer natmtot
! index to atoms and species
integer idxas(maxatoms,maxspecies)
! inverse atoms and species indices
integer idxis(maxatoms*maxspecies)
integer idxia(maxatoms*maxspecies)
! molecule is .true. is the system is an isolated molecule
logical molecule
! primcell is .true. if primitive unit cell is to be found automatically
logical primcell
! atomic positions in lattice coordinates
real(8) atposl(3,maxatoms,maxspecies)
! atomic positions in Cartesian coordinates
real(8) atposc(3,maxatoms,maxspecies)
! magnitude of random displacements added to the atomic positions
real(8) rndatposc

!----------------------------------!
!     atomic species variables     !
!----------------------------------!
! species files path
character(256) sppath
! species filenames
character(256) spfname(maxspecies)
! species name
character(256) spname(maxspecies)
! species symbol
character(256) spsymb(maxspecies)
! species nuclear charge
real(8) spzn(maxspecies)
! ptnucl is .true. if the nuclei are to be treated as point charges, if .false.
! the nuclei have a finite spherical distribution
logical ptnucl
! nuclear radius
real(8) rnucl(maxspecies)
! nuclear volume
real(8) volnucl(maxspecies)
! number of radial mesh points to nuclear radius
integer nrnucl(maxspecies)
! number of coarse radial mesh points to nuclear radius
integer nrcnucl(maxspecies)
! nuclear Coulomb potential
real(8), allocatable :: vcln(:,:)
! species electronic charge
real(8) spze(maxspecies)
! species mass
real(8) spmass(maxspecies)
! smallest radial point for each species
real(8) rminsp(maxspecies)
! effective infinity for species
real(8) rmaxsp(maxspecies)
! number of radial points to effective infinity for each species
integer nrsp(maxspecies)
! maximum nrsp over all the species
integer nrspmax
! maximum allowed states for each species
integer, parameter :: maxstsp=40
! number of states for each species
integer nstsp(maxspecies)
! maximum nstsp over all the species
integer nstspmax
! core-valence cut-off energy for species file generation
real(8) ecvcut
! semi-core-valence cut-off energy for species file generation
real(8) esccut
! state principle quantum number for each species
integer nsp(maxstsp,maxspecies)
! state l value for each species
integer lsp(maxstsp,maxspecies)
! state k value for each species
integer ksp(maxstsp,maxspecies)
! spcore is .true. if species state is core
logical spcore(maxstsp,maxspecies)
! total number of core states
integer nstcr
! state eigenvalue for each species
real(8) evalsp(maxstsp,maxspecies)
! state occupancy for each species
real(8) occsp(maxstsp,maxspecies)
! species radial mesh to effective infinity
real(8), allocatable :: rsp(:,:)
! r^l on radial mesh to muffin-tin radius
real(8), allocatable :: rlsp(:,:,:)
! species charge density
real(8), allocatable :: rhosp(:,:)
! species self-consistent potential
real(8), allocatable :: vrsp(:,:)
! exchange-correlation type for atomic species (the converged ground-state of
! the crystal does not depend on this choice)
integer xctsp(3)

end module

