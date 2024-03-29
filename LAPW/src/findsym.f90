! !INPUT/OUTPUT PARAMETERS:
!   apl1 : first set of atomic positions in lattice coordinates
!          (in,real(3,maxatoms,maxspecies))
!   apl2 : second set of atomic positions in lattice coordinates
!          (in,real(3,maxatoms,maxspecies))
!   nsym : number of symmetries (out,integer)
!   lspl : spatial rotation element in lattice point group for each symmetry
!          (out,integer(48))
!   lspn : spin rotation element in lattice point group for each symmetry
!          (out,integer(48))
!   iea  : equivalent atom index for each symmetry
!          (out,integer(iea(natmmax,nspecies,48))
! !DESCRIPTION:
!   Finds the symmetries which rotate one set of atomic positions into another.
!   Both sets of positions differ only by a translation vector and have the same
!   muffin-tin magnetic fields (stored in the global array {\tt bfcmt}). Any
!   symmetry element consists of a spatial rotation of the atomic position
!   vectors followed by a global magnetic rotation: $\{\alpha_S|\alpha_R\}$. In
!   the case of spin-orbit coupling $\alpha_S=\alpha_R$. The symmetries are
!   RETURN ed as indices of elements in the Bravais lattice point group. An
!   index to equivalent atoms is stored in the array {\tt iea}.
SUBROUTINE findsym(apl1,apl2,nsym,lspl,lspn,iea)
  USE m_atoms, ONLY: maxatoms, natoms, maxspecies, nspecies, natmmax
  USE m_spin, ONLY: bfcmt0, bfieldc0, spinpol, spinorb
  USE m_lattice, ONLY: epslat
  USE m_symmetry, ONLY: symlatd, nsymlat, symlat, symlatc, nsymlat
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(in) :: apl1(3,maxatoms,maxspecies)
  REAL(8), intent(in) :: apl2(3,maxatoms,maxspecies)
  INTEGER, intent(out) :: nsym
  INTEGER, intent(out) :: lspl(48)
  INTEGER, intent(out) :: lspn(48)
  INTEGER, intent(out) :: iea(natmmax,nspecies,48)
  ! local variables
  INTEGER :: isym,jsym,jsym0,jsym1
  INTEGER :: is,ia,ja,md
  REAL(8) :: sl(3,3),sc(3,3),v(3),t1
  ! automatic arrays
  INTEGER :: jea(natmmax,nspecies)
  REAL(8) :: apl3(3,natmmax)
  ! external functions
  REAL(8) dnrm2
  external dnrm2
  nsym=0

  ! loop over lattice symmetries (spatial rotations)
  DO isym=1,nsymlat
    ! make real copy of lattice rotation symmetry
    sl(:,:)=dble(symlat(:,:,isym))
    ! loop over species
    DO is=1,nspecies
      ! map apl1 coordinates to [0,1) and store in apl3
      DO ia=1,natoms(is)
        apl3(:,ia)=apl1(:,ia,is)
        CALL r3frac(epslat,apl3(:,ia))
      ENDDO 
      DO ja=1,natoms(is)
        ! apply lattice symmetry to atomic positions
        v(:) = sl(:,1)*apl2(1,ja,is) + sl(:,2)*apl2(2,ja,is) + sl(:,3)*apl2(3,ja,is)
        ! map coordinates to [0,1)
        CALL r3frac(epslat,v)
        ! check if atomic positions are invariant
        DO ia=1,natoms(is)
          t1=abs(apl3(1,ia)-v(1)) + abs(apl3(2,ia)-v(2)) + abs(apl3(3,ia)-v(3))
          IF(t1 < epslat) THEN 
            ! equivalent atom index
            jea(ia,is) = ja
            GOTO 10
          ENDIF 
        ENDDO 
        ! not invariant so try new spatial rotation
        GOTO 40
        10 CONTINUE
      ENDDO 
    ENDDO 
    
    ! all atomic positions invariant at this point
    jsym=1
    ! spin polarised case
    IF(spinpol) THEN 
      ! check invariance of magnetic fields under global spin rotation
      IF(spinorb) THEN 
        ! with spin-orbit coupling spin rotation equals spatial rotation
        jsym0=isym
        jsym1=isym
      ELSE 
        ! without spin-orbit coupling spin rotation independent of spatial rotation
        jsym0=1
        jsym1=nsymlat
      ENDIF 
      DO jsym=jsym0,jsym1
        ! determinant of the symmetry matrix
        md=symlatd(jsym)
        sc(:,:)=dble(md)*symlatc(:,:,jsym)
        ! rotate global field and check invariance using proper part of symmetry matrix
        v(:)=sc(:,1)*bfieldc0(1) + sc(:,2)*bfieldc0(2) + sc(:,3)*bfieldc0(3)
        t1 = abs(bfieldc0(1)-v(1)) + abs(bfieldc0(2)-v(2)) + abs(bfieldc0(3)-v(3))
        ! if not invariant try a different global spin rotation
        IF(t1.gt.epslat) goto 20
        ! rotate muffin-tin magnetic fields and check invariance
        DO is=1,nspecies
          DO ia=1,natoms(is)
            ! equivalent atom
            ja=jea(ia,is)
            v(:)=sc(:,1)*bfcmt0(1,ja,is) &
                +sc(:,2)*bfcmt0(2,ja,is) &
                +sc(:,3)*bfcmt0(3,ja,is)
            t1=abs(bfcmt0(1,ia,is)-v(1)) &
              +abs(bfcmt0(2,ia,is)-v(2)) &
              +abs(bfcmt0(3,ia,is)-v(3))
              ! if not invariant try a different global spin rotation
            IF(t1.gt.epslat) goto 20
          ENDDO 
        ENDDO 
        ! all fields invariant
        goto 30
        20 continue
      ! end loop over global spin rotations
      ENDDO 
      ! magnetic fields not invariant so try different spatial rotation
      GOTO 40
    ENDIF 
    
    30 continue
    ! everything invariant so add symmetry to set
    nsym=nsym+1
    lspl(nsym)=isym
    lspn(nsym)=jsym
    DO is=1,nspecies
      DO ia=1,natoms(is)
        iea(ia,is,nsym)=jea(ia,is)
      ENDDO 
    ENDDO 
    40 continue
    ! end loop over spatial rotations
  ENDDO 

  RETURN 
END SUBROUTINE 
