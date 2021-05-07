SUBROUTINE init_eigensystems()

  USE modmain
  IMPLICIT NONE 
  INTEGER :: lm1, lm2, lm3, l1, l2, l3, ik, jspn
  INTEGER :: m1, m2, m3
  ! external functions
  COMPLEX(8) :: gauntyry
  EXTERNAL gauntyry

  !---------------------------------------!
  !     eigenvalue equation variables     !
  !---------------------------------------!
  
  ! total number of empty states (M. Meinert)
  nempty = nint(nempty0*max(natmtot,1))
  IF( nempty < 1 ) nempty = 1
  
  ! number of first-variational states
  nstfv = int(chgval/2.d0) + nempty + 1
  
  ! overlap and Hamiltonian matrix sizes
  IF( allocated(nmat) ) deallocate(nmat)
  ALLOCATE( nmat(nspnfv,nkpt) )
  nmatmax = 0
  DO ik = 1,nkpt
    DO jspn = 1,nspnfv
      nmat(jspn,ik) = ngk(jspn,ik) + nlotot
      IF( nstfv > nmat(jspn,ik)) THEN 
        WRITE(*,*)
        WRITE(*,'("Error(init1): number of first-variational states larger than matrix size")')
        WRITE(*,'("Increase rgkmax or decrease nempty")')
        WRITE(*,*)
        STOP 
      ENDIF 
      nmatmax = max(nmatmax, nmat(jspn,ik))
    ENDDO 
  ENDDO 
  
  ! number of second-variational states
  nstsv = nstfv*nspinor
  
  ! allocate second-variational arrays
  IF( allocated(evalsv) ) DEALLOCATE(evalsv)
  ALLOCATE( evalsv(nstsv,nkpt) )
  !
  IF( allocated(occsv) ) DEALLOCATE(occsv)
  ALLOCATE( occsv(nstsv,nkpt) )
  occsv(:,:) = 0.d0
  
  ! allocate overlap and Hamiltonian integral arrays
  IF( allocated(oalo) ) DEALLOCATE(oalo)
  ALLOCATE( oalo(apwordmax,nlomax,natmtot) )
  !
  IF( allocated(ololo) ) DEALLOCATE(ololo)
  ALLOCATE(ololo(nlomax,nlomax,natmtot))
  IF( allocated(haa) ) DEALLOCATE(haa)
  ALLOCATE( haa(lmmaxo,apwordmax,0:lmaxapw,apwordmax,0:lmaxapw,natmtot) )
  !
  IF( allocated(hloa) ) DEALLOCATE(hloa)
  allocate(hloa(lmmaxo,apwordmax,0:lmaxapw,nlomax,natmtot))
  !
  IF( allocated(hlolo)) DEALLOCATE(hlolo)
  ALLOCATE( hlolo(lmmaxo,nlomax,nlomax,natmtot) )
  
  ! allocate and generate complex Gaunt coefficient array
  IF( allocated(gntyry) ) deallocate(gntyry)
  ALLOCATE( gntyry(lmmaxo,lmmaxapw,lmmaxapw) )
  DO l1 = 0,lmaxapw
    DO m1 = -l1,l1
      lm1 = idxlm(l1,m1)
      DO l3 = 0,lmaxapw
        DO m3=-l3,l3
          lm3 = idxlm(l3,m3)
          DO l2 = 0,lmaxo
            DO m2 = -l2,l2
              lm2 = idxlm(l2,m2)
              gntyry(lm2,lm3,lm1) = gauntyry(l1,l2,l3,m1,m2,m3)
            ENDDO 
          ENDDO 
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO 
  ! write to VARIABLES.OUT
  !call writevars('nempty',iv=nempty)
  !call writevars('nstfv',iv=nstfv)
  !call writevars('nlotot',iv=nlotot)
  !call writevars('nstsv',iv=nstsv)

  RETURN 

END SUBROUTINE 
