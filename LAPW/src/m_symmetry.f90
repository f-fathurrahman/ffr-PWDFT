module m_symmetry

!----------------------------!
!     symmetry variables     !
!----------------------------!

! type of symmetry allowed for the crystal
!  0 : only the identity element is used
!  1 : full symmetry group is used
!  2 : only symmorphic symmetries are allowed
integer symtype

! number of Bravais lattice point group symmetries
integer nsymlat

! Bravais lattice point group symmetries
integer symlat(3,3,48)

! determinants of lattice symmetry matrices (1 or -1)
integer symlatd(48)

! index to inverses of the lattice symmetries
integer isymlat(48)

! lattice point group symmetries in Cartesian coordinates
real(8) symlatc(3,3,48)

! tshift is .true. if atomic basis is allowed to be shifted
logical tshift

! tsyminv is .true. if the crystal has inversion symmetry
logical tsyminv

! maximum of symmetries allowed
integer, parameter :: maxsymcrys=192

! number of crystal symmetries
integer nsymcrys

! crystal symmetry translation vector in lattice and Cartesian coordinates
real(8) vtlsymc(3,maxsymcrys)
real(8) vtcsymc(3,maxsymcrys)
! tv0symc is .true. if the translation vector is zero
logical tv0symc(maxsymcrys)
! spatial rotation element in lattice point group for each crystal symmetry
integer lsplsymc(maxsymcrys)
! global spin rotation element in lattice point group for each crystal symmetry
integer lspnsymc(maxsymcrys)
! equivalent atom index for each crystal symmetry
integer, allocatable :: ieqatom(:,:,:)
! eqatoms(ia,ja,is) is .true. if atoms ia and ja are equivalent
logical, allocatable :: eqatoms(:,:,:)
! number of site symmetries
integer, allocatable :: nsymsite(:)
! site symmetry spatial rotation element in lattice point group
integer, allocatable :: lsplsyms(:,:)
! site symmetry global spin rotation element in lattice point group
integer, allocatable :: lspnsyms(:,:)

end module 