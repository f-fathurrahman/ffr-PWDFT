SUBROUTINE init_PsPot()

  USE m_atoms, ONLY : SpeciesSymbols, Nspecies, AtomicValences, &
                      atm2species, Natoms
  USE m_PsPot, ONLY : Ps_HGH_Params, PsPot_FilePath, w_NL, PsPot_Dir, &
                      NbetaNL
  USE m_Ps_HGH, ONLY : init_Ps_HGH_Params

  IMPLICIT NONE 
  INTEGER :: isp, iprj, l, ia, m, ibeta

  ! Use default
  ALLOCATE( PsPot_FilePath(Nspecies) )
  ALLOCATE( Ps_HGH_Params(Nspecies) )

  DO isp = 1,Nspecies

    PsPot_FilePath(isp) = trim(PsPot_Dir) // trim(SpeciesSymbols(isp)) // '.hgh'
    CALL init_Ps_HGH_Params( Ps_HGH_Params(isp), PsPot_FilePath(isp) )
    
    AtomicValences(isp) = Ps_HGH_Params(isp)%zval

  ENDDO 

  DO ia = 1,Natoms
    isp = atm2species(ia)
    DO l = 0,Ps_HGH_Params(isp)%lmax
      DO iprj = 1,Ps_HGH_Params(isp)%Nproj_l(l)
        DO m = -l,l
          NbetaNL = NbetaNL + 1
        ENDDO ! m
      ENDDO ! iprj
    ENDDO ! l
  ENDDO 
  WRITE(*,*) 'NbetaNL = ', NbetaNL

  ALLOCATE( w_NL(NbetaNL) )

  ibeta = 0
  DO ia = 1,Natoms
    isp = atm2species(ia)
    DO l = 0,Ps_HGH_Params(isp)%lmax
      DO iprj = 1,Ps_HGH_Params(isp)%Nproj_l(l)
        DO m = -l,l
          ibeta = ibeta + 1
          w_NL(ibeta) = Ps_HGH_Params(isp)%h(l,iprj,iprj)
        ENDDO ! m
      ENDDO ! iprj
    ENDDO ! l
  ENDDO 
  WRITE(*,*) 'w_NL = ', w_NL

END SUBROUTINE 

