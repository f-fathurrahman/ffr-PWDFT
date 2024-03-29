!---------------------------------
SUBROUTINE init_idx_atom_species()
!---------------------------------
  USE m_atoms, ONLY: nspecies, idxas, idxis, idxia, natmmax, natmtot, natoms
  IMPLICIT NONE 
  INTEGER :: is, ia, ias

  WRITE(*,*) 'Setting up index to atoms and species'

  !------------------------------------!
  !     index to atoms and species     !
  !------------------------------------!
  natmmax = 0
  ias = 0
  DO is = 1,nspecies
    DO ia = 1,natoms(is)
      ias = ias + 1
      idxas(ia,is) = ias
      idxis(ias) = is
      idxia(ias) = ia
    ENDDO 
    ! maximum number of atoms over all species
    natmmax = max(natmmax,natoms(is))
  ENDDO 
  
  ! total number of atoms
  natmtot = ias
  
  ! number of phonon branches
  !nbph = 3*natmtot
  
  ! write to VARIABLES.OUT
  !call writevars('nspecies',iv=nspecies)
  !call writevars('natoms',nv=nspecies,iva=natoms)
  !call writevars('spsymb',nv=nspecies,sva=spsymb)
  !call writevars('spname',nv=nspecies,sva=spname)
  !call writevars('spzn',nv=nspecies,rva=spzn)

END SUBROUTINE 
