MODULE m_atomic_species

USE m_atoms, ONLY: maxspecies

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
REAL(8) spzn(maxspecies)
! ptnucl is .true. if the nuclei are to be treated as point charges, if .false.
! the nuclei have a finite spherical distribution
logical ptnucl
! nuclear radius
REAL(8) rnucl(maxspecies)
! nuclear volume
REAL(8) volnucl(maxspecies)
! number of radial mesh points to nuclear radius
INTEGER nrnucl(maxspecies)
! number of coarse radial mesh points to nuclear radius
INTEGER nrcnucl(maxspecies)
! nuclear Coulomb potential
REAL(8), ALLOCATABLE :: vcln(:,:)
! species electronic charge
REAL(8) spze(maxspecies)
! species mass
REAL(8) spmass(maxspecies)
! smallest radial point for each species
REAL(8) rminsp(maxspecies)
! effective infinity for species
REAL(8) rmaxsp(maxspecies)
! number of radial points to effective infinity for each species
INTEGER nrsp(maxspecies)
! maximum nrsp over all the species
INTEGER nrspmax
! maximum allowed states for each species
INTEGER, parameter :: maxstsp=40
! number of states for each species
INTEGER nstsp(maxspecies)
! maximum nstsp over all the species
INTEGER nstspmax
! core-valence cut-off energy for species file generation
REAL(8) ecvcut
! semi-core-valence cut-off energy for species file generation
REAL(8) esccut
! state principle quantum number for each species
INTEGER nsp(maxstsp,maxspecies)
! state l value for each species
INTEGER lsp(maxstsp,maxspecies)
! state k value for each species
INTEGER ksp(maxstsp,maxspecies)
! spcore is .true. if species state is core
LOGICAL spcore(maxstsp,maxspecies)
! total number of core states
INTEGER nstcr
! state eigenvalue for each species
REAL(8) evalsp(maxstsp,maxspecies)
! state occupancy for each species
REAL(8) occsp(maxstsp,maxspecies)
! species radial mesh to effective infinity
REAL(8), ALLOCATABLE :: rsp(:,:)
! r^l on radial mesh to muffin-tin radius
REAL(8), ALLOCATABLE :: rlsp(:,:,:)
! species charge density
REAL(8), ALLOCATABLE :: rhosp(:,:)
! species self-consistent potential
REAL(8), ALLOCATABLE :: vrsp(:,:)
! exchange-correlation type for atomic species (the converged ground-state of
! the crystal does not depend on this choice)
INTEGER xctsp(3)

END MODULE 