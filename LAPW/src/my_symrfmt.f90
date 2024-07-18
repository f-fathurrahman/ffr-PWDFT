SUBROUTINE my_symrfmt(nr, nri, np, ld, rfmt)
  !
  ! The input/output rfmt is for all atoms (2d array with size npmtmax x natmtom)
  !
  USE m_atoms, ONLY: nspecies, natmtot, natoms, idxas, natmmax
  USE m_symmetry, ONLY: nsymcrys, lsplsymc, ieqatom, isymlat, symlatc
  !
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr(nspecies), nri(nspecies), np(nspecies)
  INTEGER, intent(in) :: ld
  REAL(8), intent(inout) :: rfmt(ld,natmtot)
  !
  ! local variables
  INTEGER is,ia,ja,ias,jas
  INTEGER isym,lspl
  REAL(8) t0
  ! automatic arrays
  LOGICAL done(natmmax)
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rfmt1(:,:), rfmt2(:)
  !
  ALLOCATE( rfmt1(ld,natmmax), rfmt2(ld) )
  t0 = 1.d0/dble(nsymcrys)
  !
  DO is = 1,nspecies
    ! make a copy of the input function
    DO ia = 1,natoms(is)
      ias = idxas(ia,is)
      CALL dcopy(np(is), rfmt(:,ias),1, rfmt1(:,ia), 1)
    ENDDO 
    done(:) = .false.
    !
    ! loop over atoms
    DO ia = 1,natoms(is)
      !
      IF( done(ia) ) cycle
      !
      ias = idxas(ia,is)
      rfmt(1:np(is),ias) = 0.d0
      ! loop over crystal symmetries
      DO isym = 1,nsymcrys
        !
        ! index to spatial rotation lattice symmetry
        lspl = lsplsymc(isym)
        !
        ! equivalent atom index (symmetry rotates atom ja into atom ia)
        ja = ieqatom(ia,is,isym)
        !
        ! apply the rotation to the muffin-tin function
        CALL my_rotrfmt(symlatc(:,:,lspl), nr(is), nri(is), rfmt1(:,ja), rfmt2)
        !
        ! accumulate in original function array
        rfmt(1:np(is),ias) = rfmt(1:np(is),ias) + rfmt2(1:np(is))
      ENDDO 
      ! normalise
      CALL dscal(np(is), t0, rfmt(:,ias), 1)
      done(ia) = .true.
      ! rotate into equivalent atoms
      DO isym = 1,nsymcrys
        ja = ieqatom(ia, is, isym)
        IF( done(ja) ) cycle
        jas = idxas(ja,is)
        !
        ! inverse symmetry (which rotates atom ia into atom ja)
        lspl = isymlat( lsplsymc(isym) )
        !
        ! rotate symmetrised function into equivalent muffin-tin
        CALL my_rotrfmt( symlatc(:,:,lspl), nr(is), nri(is), rfmt(:,ias), rfmt(:,jas) )
        done(ja) = .true.
      ENDDO 
    ENDDO 
  ENDDO 
  DEALLOCATE(rfmt1,rfmt2)
  RETURN 
END SUBROUTINE 

