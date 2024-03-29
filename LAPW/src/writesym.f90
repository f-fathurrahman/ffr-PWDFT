SUBROUTINE writesym
  USE modmain
! !DESCRIPTION:
!   Outputs the Bravais, crystal and site symmetry matrices to files
!   {\tt SYMLAT.OUT}, {\tt SYMCRYS.OUT} and {\tt SYMSITE.OUT}, respectively.
!   Also writes out equivalent atoms and related crystal symmetries to
!   {\tt EQATOMS.OUT}.
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ia,ja,ias,i
  INTEGER :: isym,lspl,lspn

  ! output the Bravais lattice symmetries
  open(50,file='SYMLAT'//trim(filext),form='FORMATTED')
  WRITE(50,'(I4," : nsymlat")') nsymlat
  DO isym=1,nsymlat
    WRITE(50,*)
    WRITE(50,'(I4)') isym
    DO i=1,3
      WRITE(50,'(3I4)') symlat(i,:,isym)
    ENDDO 
  ENDDO 
  close(50)

  ! output the crystal symmetries
  open(50,file='SYMCRYS'//trim(filext),form='FORMATTED')
  WRITE(50,*)
  WRITE(50,'("(translation vectors and rotation matrices are in lattice &
   &coordinates)")')
  WRITE(50,*)
  WRITE(50,'(I4," : nsymcrys")') nsymcrys
  DO isym=1,nsymcrys
    WRITE(50,*)
    WRITE(50,'("Crystal symmetry : ",I4)') isym
    WRITE(50,'(" spatial translation :")')
    WRITE(50,'(3G18.10)') vtlsymc(:,isym)
    WRITE(50,'(" spatial rotation :")')
    lspl=lsplsymc(isym)
    DO i=1,3
      WRITE(50,'(3I4)') symlat(i,:,lspl)
    ENDDO 
    WRITE(50,'(" global spin rotation :")')
    lspn=lspnsymc(isym)
    DO i=1,3
      WRITE(50,'(3I4)') symlat(i,:,lspn)
    ENDDO 
  ENDDO 
  close(50)

  ! output the site symmetries
  open(50,file='SYMSITE'//trim(filext),form='FORMATTED')
  WRITE(50,*)
  WRITE(50,'("(rotation matrices are in lattice coordinates)")')
  DO is=1,nspecies
    DO ia=1,natoms(is)
      ias=idxas(ia,is)
      WRITE(50,*)
      WRITE(50,*)
      WRITE(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
      WRITE(50,'(I4," : nsymsite")') nsymsite(ias)
      DO isym=1,nsymsite(ias)
        WRITE(50,*)
        WRITE(50,'(" Site symmetry : ",I4)') isym
        WRITE(50,'("  spatial rotation :")')
        lspl=lsplsyms(isym,ias)
        DO i=1,3
          WRITE(50,'(3I4)') symlat(i,:,lspl)
        ENDDO 
        WRITE(50,'("  global spin rotation :")')
        lspn=lspnsyms(isym,ias)
        DO i=1,3
          WRITE(50,'(3I4)') symlat(i,:,lspn)
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO 
  close(50)

  ! output the equivalent atoms and related symmetries
  open(50,file='EQATOMS'//trim(filext),form='FORMATTED')
  DO is=1,nspecies
    WRITE(50,*)
    WRITE(50,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
    DO ia=1,natoms(is)
      WRITE(50,'(" atom ",I4," is equivalent to atom(s)")') ia
      i=0
      DO ja=1,natoms(is)
        IF(eqatoms(ia,ja,is)) THEN 
          IF((i.gt.0).and.(mod(i,20).eq.0)) WRITE(50,*)
          WRITE(50,'(I4)',advance='NO') ja
          i=i+1
        ENDIF 
      ENDDO 
      WRITE(50,*)
    ENDDO 
  ENDDO 
  close(50)
  
  RETURN 
END SUBROUTINE 
