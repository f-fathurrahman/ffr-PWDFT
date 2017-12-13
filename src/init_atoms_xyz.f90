

SUBROUTINE init_atoms_xyz( fil_xyz )

  USE m_constants, ONLY: ANG2BOHR
  USE m_atoms

  IMPLICIT NONE
  ! Argument
  CHARACTER(*) :: fil_xyz

  ! Local
  INTEGER, PARAMETER :: unitxyz=123
  INTEGER :: ios, ia, k1, k2, isp, idx1
  CHARACTER(5), ALLOCATABLE :: atmSymb(:)

  !
  OPEN( UNIT=unitxyz, FILE=fil_xyz, ACTION='read', STATUS='old', &
    FORM='formatted', IOSTAT=ios )
  IF( ios /= 0 ) THEN 
    WRITE(*,*) 'ERROR reading file: ', trim(fil_xyz)
    STOP 
  ENDIF 

  !
  ! Read
  !
  READ( unitxyz, * ) Natoms
  READ( unitxyz, * )

  ALLOCATE( atmSymb( Natoms ) )

  ALLOCATE( AtomicCoords( 3, Natoms ) )

  DO ia=1,natoms
    READ( unitxyz, * ) atmSymb(ia), &
                       AtomicCoords(1,ia), AtomicCoords(2,ia), AtomicCoords(3,ia)
    !WRITE(*,'(1x,A,3G18.10)') trim(atoms%atmSymb(ia)), atoms%positions(1,ia), &
    !  atoms%positions(2,ia), atoms%positions(3,ia)
  ENDDO
  CLOSE(unitxyz)

  !
  ! convert from angstrom to bohr
  AtomicCoords(:,:) = AtomicCoords(:,:) * ANG2BOHR

  ! Determine number of species
  Nspecies = 0
  DO ia=1,Natoms
    k2 = 0
    DO k1=1,ia-1
      IF( atmSymb(k1) == atmSymb(ia) ) k2=1
    ENDDO
    ! find different
    IF(k2==0) THEN
      Nspecies = Nspecies + 1
    ENDIF
  ENDDO

  ALLOCATE( SpeciesSymbols(NSPECIES) )

  idx1 = 0
  DO ia=1,Natoms
    k2 = 0
    DO k1=1,ia-1
      IF( atmSymb(k1) == atmSymb(ia) ) k2=1
    ENDDO
    ! Found different species
    IF(k2==0) THEN
      idx1 = idx1 + 1
      SpeciesSymbols(idx1) = atmSymb(ia)
    ENDIF
  ENDDO
  
  ! Mapping of atoms to species index
  ALLOCATE( atm2species(NATOMS) )
  DO ia=1,Natoms
    DO isp=1,Nspecies
      IF( atmSymb(ia) == SpeciesSymbols(isp) ) THEN
        atm2Species(ia) = isp
      ENDIF
    ENDDO 
  ENDDO

  ! 
  ALLOCATE( AtomicValences(Nspecies) )
  AtomicValences(:) = 0.d0  ! NOTE: They should be set by pseudopotentials or manually

END SUBROUTINE

