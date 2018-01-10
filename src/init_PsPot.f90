SUBROUTINE init_PsPot()

  USE m_atoms, ONLY : SpeciesSymbols, Nspecies, AtomicValences, &
                      atm2species, Natoms
  USE m_PsPot, ONLY : Ps_HGH_Params, PsPot_FilePath, w_NL, PsPot_Dir, &
                      NbetaNL, PsPot_lll, PsPot_ipr, NprojTotMax, NprojTot => SpecNprojTot
  USE m_Ps_HGH, ONLY : init_Ps_HGH_Params

  IMPLICIT NONE 
  INTEGER :: isp, iprj, l, ia, m, ibeta
  INTEGER :: iv, ii

  ! Use default
  ALLOCATE( PsPot_FilePath(Nspecies) )
  ALLOCATE( Ps_HGH_Params(Nspecies) )

  DO isp = 1,Nspecies

    PsPot_FilePath(isp) = trim(PsPot_Dir) // trim(SpeciesSymbols(isp)) // '.gth'

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
  IF( NbetaNL > 0 ) THEN 
    WRITE(*,*) 'w_NL = ', w_NL
  ENDIF 

  ALLOCATE( NprojTot(Nspecies) )
  NprojTotMax = 0
  DO isp = 1,Nspecies
    NprojTot(isp) = sum( Ps_HGH_Params(isp)%Nproj_l(:) )
    NprojTotMax = max( NprojTotMax, NprojTot(isp) )
    !WRITE(*,*) 'isp, NprojTot = ', isp, NprojTot(isp)
  ENDDO
  !WRITE(*,*) 'NprojTotMax = ', NprojTotMax


  ALLOCATE( PsPot_lll(Nspecies,NprojTotMax) )
  ALLOCATE( PsPot_ipr(Nspecies,NprojTotMax) )

  PsPot_lll(:,:) = -1
  PsPot_ipr(:,:) = -1

  ! For mk_ffnl_gth
  DO isp = 1,Nspecies
    IF( Ps_HGH_Params(isp)%lmax /= -1 ) THEN 
      iv = 0
      DO l = 0,Ps_HGH_Params(isp)%lmax
        DO ii = 1,Ps_HGH_Params(isp)%Nproj_l(l)
          iv = iv + 1
          PsPot_lll(isp,iv) = l
          PsPot_ipr(isp,iv) = ii
        ENDDO 
      ENDDO 
    ENDIF 
!    WRITE(*,*) 'isp, iv = ', isp, iv
  ENDDO 

!  DO isp = 1,Nspecies
!    IF( NprojTot(isp) == 0 ) THEN 
!      WRITE(*,*) 'Species ', isp, ' does not have NL PS projectors.'
!    ENDIF 
!    DO iprj = 1,NprojTot(isp)
!      WRITE(*,*) iprj, PsPot_lll(isp,iprj), PsPot_ipr(isp,iprj)
!    ENDDO 
!  ENDDO 


END SUBROUTINE 

