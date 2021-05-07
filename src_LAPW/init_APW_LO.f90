SUBROUTINE init_APW_LO()

  USE modmain
  IMPLICIT NONE 
  INTEGER :: ia, ias, ilo, l1, is, io

  !---------------------------------!
  !     APWs and local-orbitals     !
  !---------------------------------!
  apwordmax = 0
  lorbordmax = 0
  nlomax = 0
  lolmax = 0
  DO is = 1,nspecies
    lmoapw(is) = 0
    DO l1 = 0,lmaxapw
      ! find the maximum APW order
      apwordmax = max(apwordmax,apword(l1,is))
      ! find total number of APW coefficients (l, m and order)
      lmoapw(is) = lmoapw(is)+(2*l1+1)*apword(l1,is)
    ENDDO 
    ! find the maximum number of local-orbitals
    nlomax = max(nlomax, nlorb(is))
    ! find the maximum local-orbital order and angular momentum
    DO ilo = 1,nlorb(is)
      lolmax = max( lolmax, lorbl(ilo,is) )
      lorbordmax = max( lorbordmax, lorbord(ilo,is) )
    ENDDO 
  ENDDO 
  lolmmax = (lolmax + 1)**2
  
  ! polynomial order used for APW and local-orbital radial derivatives
  npapw = max(apwordmax+1, 4)
  nplorb = max(lorbordmax+1, 4)
  
  ! set the APW and local-orbital linearisation energies to the default
  IF( allocated(apwe) ) DEALLOCATE(apwe)
  allocate( apwe(apwordmax, 0:lmaxapw,natmtot) )
  !
  IF( allocated(lorbe) ) DEALLOCATE(lorbe)
  ALLOCATE( lorbe(lorbordmax,maxlorb,natmtot) )
  !
  DO is = 1,nspecies
    DO l1 = 0,lmaxapw
      DO io = 1,apword(l1,is)
        do ia = 1,natoms(is)
          ias = idxas(ia,is)
          apwe(io,l1,ias) = apwe0(io,l1,is)
        ENDDO 
      ENDDO 
    ENDDO 
    DO ilo = 1,nlorb(is)
      DO io = 1,lorbord(ilo,is)
        DO ia = 1,natoms(is)
          ias=idxas(ia,is)
          lorbe(io,ilo,ias) = lorbe0(io,ilo,is)
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO 
  
  ! generate the local-orbital index
  CALL genidxlo()
  
  ! allocate radial function arrays
  IF( allocated(apwfr) ) DEALLOCATE(apwfr)
  allocate(apwfr(nrmtmax,2,apwordmax,0:lmaxapw,natmtot))
  !
  IF( allocated(apwdfr) ) DEALLOCATE(apwdfr)
  ALLOCATE(apwdfr(apwordmax,0:lmaxapw,natmtot))
  !
  IF( allocated(lofr) ) DEALLOCATE(lofr)
  ALLOCATE( lofr(nrmtmax,2,nlomax,natmtot) )



END SUBROUTINE 

