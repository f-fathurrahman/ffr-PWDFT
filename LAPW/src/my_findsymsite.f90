SUBROUTINE my_findsymsite()

  use m_atoms, only: natoms, natmtot, nspecies, maxatoms, maxspecies, idxas, natmmax, atposl
  use m_symmetry, only: lspnsyms, nsymsite, lsplsyms

  IMPLICIT NONE 

  ! local variables
  INTEGER is,js,ia,ja,ias
  REAL(8) apl(3,maxatoms,maxspecies)
  ! automatic arrays
  REAL(8) iea(natmmax,nspecies,48)
  
  write(*,*) '--------------------'
  write(*,*) 'Enter my_findsymsite'
  write(*,*) '--------------------'

  ! allocate the site symmetry arrays
  IF(allocated(nsymsite)) DEALLOCATE(nsymsite)
  ALLOCATE(nsymsite(natmtot))
  !
  IF(allocated(lsplsyms)) DEALLOCATE(lsplsyms)
  ALLOCATE(lsplsyms(48,natmtot))
  !
  IF(allocated(lspnsyms)) DEALLOCATE(lspnsyms)
  ALLOCATE(lspnsyms(48,natmtot))
  
  DO is = 1,nspecies
    DO ia = 1,natoms(is)
      ias = idxas(ia,is)
      DO js = 1,nspecies
        DO ja = 1,natoms(js)
          apl(:,ja,js) = atposl(:,ja,js) - atposl(:,ia,is)
        ENDDO 
      ENDDO 
      ! apl contains distance between atoms ia and ja
      !
      ! why the first and second arguments are the same?
      CALL my_findsym(apl, apl, nsymsite(ias), lsplsyms(:,ias), lspnsyms(:,ias), iea)
      write(*,'(1x,A,2I4)') 'ias, nsymsite(ias) = ', ias, nsymsite(ias)
    ENDDO 
  ENDDO 

  write(*,*) '--------------------'
  write(*,*) 'Exit my_findsymsite'
  write(*,*) '--------------------'

  RETURN 

END SUBROUTINE 

