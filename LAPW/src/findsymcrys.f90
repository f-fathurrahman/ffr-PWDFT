! !DESCRIPTION:
!   Finds the complete set of symmetries which leave the crystal structure
!   (including the magnetic fields) invariant. A crystal symmetry is of the form
!   $\{\alpha_S|\alpha_R|{\bf t}\}$, where ${\bf t}$ is a translation vector,
!   $\alpha_R$ is a spatial rotation operation and $\alpha_S$ is a global spin
!   rotation. Note that the order of operations is important and defined to be
!   from right to left, i.e. translation followed by spatial rotation followed
!   by spin rotation. In the case of spin-orbit coupling $\alpha_S=\alpha_R$. In
!   order to determine the translation vectors, the entire atomic basis is
!   shifted so that the first atom in the smallest set of atoms of the same
!   species is at the origin. Then all displacement vectors between atoms in
!   this set are checked as possible symmetry translations. If the global
!   variable {\tt tshift} is set to {\tt .false.} THEN  the shift is not
!   performed. See L. M. Sandratskii and P. G. Guletskii, {\it J. Phys. F: Met.
!   Phys.} {\bf 16}, L43 (1986) and the routine {\tt findsym}.
SUBROUTINE findsymcrys()
  ! !USES:
  USE m_atoms, ONLY: nspecies, atposl, maxspecies, maxatoms, atposc, natoms, &
               natmtot, natmmax, nspecies
  USE m_lattice, ONLY: epslat, avec
  USE m_electric_vector_pot, ONLY: efieldc, tefield
  USE m_symmetry, only: tv0symc, tsyminv, vtlsymc, lsplsymc, ieqatom, &
                     lspnsymc, isymlat, eqatoms, symtype, &
                     nsymcrys, tshift, maxsymcrys, &
                     vtcsymc, symlat

  IMPLICIT NONE 
  ! local variables
  INTEGER :: ia,ja,is,js
  INTEGER :: isym,nsym,i,n
  INTEGER :: lspl(48),lspn(48),ilspl
  REAL(8) :: v0(3),v1(3),v2(3),t1
  REAL(8) :: apl(3,maxatoms,maxspecies)
  ! ALLOCATABLE arrays
  INTEGER, ALLOCATABLE :: iea(:,:,:)
  REAL(8), ALLOCATABLE :: vtl(:,:)

  ! allocate local array
  ALLOCATE(iea(natmmax,nspecies,48))

  ! allocate equivalent atom arrays
  IF(allocated(ieqatom)) DEALLOCATE(ieqatom)
  ALLOCATE(ieqatom(natmmax,nspecies,maxsymcrys))
  IF(allocated(eqatoms)) DEALLOCATE(eqatoms)
  ALLOCATE(eqatoms(natmmax,natmmax,nspecies))

  ! store position of first atom
  IF(natmtot > 0) v0(:)=atposl(:,1,1)

  ! find the smallest set of atoms
  is=1
  DO js=1,nspecies
    IF(natoms(js).lt.natoms(is)) is=js
  ENDDO 

  IF((tshift).and.(natmtot > 0)) THEN 
    ! shift basis so that the first atom in the smallest atom set is at the origin
    v1(:)=atposl(:,1,is)
    DO js=1,nspecies
      DO ia=1,natoms(js)
        ! shift atom
        atposl(:,ia,js)=atposl(:,ia,js)-v1(:)
        ! map lattice coordinates back to [0,1)
        CALL r3frac(epslat,atposl(:,ia,js))
        ! determine the new Cartesian coordinates
        CALL r3mv(avec,atposl(:,ia,js),atposc(:,ia,js))
      ENDDO 
    ENDDO 
  ENDIF 

  ! determine possible translation vectors from smallest set of atoms
  n=max(natoms(is)*natoms(is),1)
  ALLOCATE(vtl(3,n))
  n=1
  vtl(:,1)=0.d0
  DO ia=1,natoms(is)
    DO ja=2,natoms(is)
      ! compute difference between two atom vectors
      v1(:)=atposl(:,ia,is)-atposl(:,ja,is)
      ! map lattice coordinates to [0,1)
      CALL r3frac(epslat,v1)
      ! check if vector has any component along electric field
      IF(tefield) THEN 
        CALL r3mv(avec,v1,v2)
        t1=efieldc(1)*v2(1)+efieldc(2)*v2(2)+efieldc(3)*v2(3)
        IF(abs(t1) > epslat) goto 10
      ENDIF 
      DO i=1,n
        t1=abs(vtl(1,i)-v1(1))+abs(vtl(2,i)-v1(2))+abs(vtl(3,i)-v1(3))
        IF(t1.lt.epslat) goto 10
      ENDDO 
      n=n+1
      vtl(:,n)=v1(:)
      10 CONTINUE
    ENDDO 
  ENDDO 

  ! no translations required when symtype=0,2 (F. Cricchio)
  IF(symtype /= 1) n=1
  eqatoms(:,:,:)=.false.
  nsymcrys=0

  ! loop over all possible translations
  DO i=1,n
    ! construct new array with translated positions
    DO is=1,nspecies
      DO ia=1,natoms(is)
        apl(:,ia,is)=atposl(:,ia,is)+vtl(:,i)
      ENDDO 
    ENDDO 
    ! find the symmetries for current translation
    CALL findsym(atposl,apl,nsym,lspl,lspn,iea)
    DO isym = 1,nsym
      nsymcrys = nsymcrys + 1
      IF(nsymcrys > maxsymcrys) THEN 
        WRITE(*,*)
        WRITE(*,'("Error(findsymcrys): too many crystal symmetries")')
        WRITE(*,'(" Adjust maxsymcrys in modmain and recompile code")')
        WRITE(*,*)
        stop
      ENDIF 
      vtlsymc(:,nsymcrys)=vtl(:,i)
      lsplsymc(nsymcrys)=lspl(isym)
      lspnsymc(nsymcrys)=lspn(isym)
      DO is=1,nspecies
        DO ia=1,natoms(is)
          ja=iea(ia,is,isym)
          ieqatom(ia,is,nsymcrys)=ja
          eqatoms(ia,ja,is)=.true.
          eqatoms(ja,ia,is)=.true.
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO 

tsyminv=.false.

DO isym=1,nsymcrys
! check if inversion symmetry is present
  i=lsplsymc(isym)
  IF(all(symlat(:,:,i) == -symlat(:,:,1))) THEN 
    tsyminv=.true.
! make inversion the second symmetry element (the identity is the first)
    v1(:)=vtlsymc(:,isym); vtlsymc(:,isym)=vtlsymc(:,2); vtlsymc(:,2)=v1(:)
    i=lsplsymc(isym); lsplsymc(isym)=lsplsymc(2); lsplsymc(2)=i
    i=lspnsymc(isym); lspnsymc(isym)=lspnsymc(2); lspnsymc(2)=i
    DO is=1,nspecies
      DO ia=1,natoms(is)
        i=ieqatom(ia,is,isym)
        ieqatom(ia,is,isym)=ieqatom(ia,is,2)
        ieqatom(ia,is,2)=i
      ENDDO 
    ENDDO 
    goto 20
  ENDIF 
ENDDO 
20 continue

! if inversion exists THEN  shift basis so that inversion center is at origin
IF(tsyminv.and.tshift) THEN 
  v1(:)=v1(:)/2.d0
  DO is=1,nspecies
    DO ia=1,natoms(is)
      ! shift atom
      atposl(:,ia,is)=atposl(:,ia,is)+v1(:)
      ! map lattice coordinates back to [0,1)
      CALL r3frac(epslat,atposl(:,ia,is))
      ! map lattice coordinates to [-0.5,0.5)
      DO i=1,3
        IF(atposl(i,ia,is) > 0.5d0) atposl(i,ia,is)=atposl(i,ia,is)-1.d0
      ENDDO 
      ! determine the new Cartesian coordinates
      CALL r3mv(avec,atposl(:,ia,is),atposc(:,ia,is))
    ENDDO 
  ENDDO 
  ! recalculate crystal symmetry translation vectors
  DO isym=1,nsymcrys
    ilspl=isymlat(lsplsymc(isym))
    v2(:)=symlat(:,1,ilspl)*v1(1) &
         +symlat(:,2,ilspl)*v1(2) &
         +symlat(:,3,ilspl)*v1(3)
    vtlsymc(:,isym)=vtlsymc(:,isym)-v1(:)+v2(:)
    CALL r3frac(epslat,vtlsymc(:,isym))
  ENDDO 
ENDIF 

! translation vector in Cartesian coordinates
DO isym=1,nsymcrys
  CALL r3mv(avec,vtlsymc(:,isym),vtcsymc(:,isym))
ENDDO 

! set flag for zero translation vector
DO isym=1,nsymcrys
  t1=abs(vtlsymc(1,isym))+abs(vtlsymc(2,isym))+abs(vtlsymc(3,isym))
  IF(t1.lt.epslat) THEN 
    tv0symc(isym)=.true.
  else
    tv0symc(isym)=.false.
  ENDIF 
ENDDO 

! check inversion does not include a translation
IF(tsyminv) THEN 
  IF(.not.tv0symc(2)) tsyminv=.false.
ENDIF 

IF(natmtot > 0) THEN 
  v1(:)=atposl(:,1,1)-v0(:)
  t1=abs(v1(1))+abs(v1(2))+abs(v1(3))
  if( t1 > epslat ) THEN 
    WRITE(*,*)
    WRITE(*,'("Info(findsymcrys): atomic basis shift (lattice) :")')
    WRITE(*,'(3G18.10)') v1(:)
    WRITE(*,'("See GEOMETRY.OUT for new atomic positions")')
  ENDIF 
ENDIF 

DEALLOCATE(iea,vtl)
RETURN 

END SUBROUTINE 
