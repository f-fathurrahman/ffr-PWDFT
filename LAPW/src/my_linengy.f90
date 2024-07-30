!----------------------
SUBROUTINE my_linengy()
!----------------------
  USE m_atoms, ONLY: natmmax, idxas, natoms, nspecies
  USE m_apwlo, ONLY: nlorb, lorbord, lorbe, apword, apwe, &
               lorbl, lorbe0, lorbe, lorbve, apwe0, apwve, dlefe, &
               autolinengy, demaxbnd, epsband
  USE m_symmetry, ONLY: eqatoms
  USE m_muffin_tins, ONLY: nrmtmax, nrmti, nrmt, lmmaxi, lmmaxo, lmaxapw, rlmt
  USE m_constants, ONLY: y00, solsc
  USE m_density_pot_xc, ONLY: vsmt
  USE m_convergence, ONLY: iscl
  USE m_states, ONLY: efermi
  !
  IMPLICIT NONE 
  ! local variables
  LOGICAL :: fnd
  INTEGER :: is,ia,ja,ias,jas
  INTEGER :: nr,nri,iro,ir,l,i
  INTEGER :: ilo,io,jo,nnf
  ! automatic arrays
  LOGICAL :: done(natmmax)
  REAL(8) :: vr(nrmtmax)
  
  nnf = 0
  ! begin loops over atoms and species
  DO is = 1,nspecies
    nr = nrmt(is)
    nri = nrmti(is)
    iro = nri + 1
    done(:) = .false.
    DO ia = 1,natoms(is)
      IF(done(ia)) CYCLE
      ias = idxas(ia,is)
      i = 1
      DO ir = 1,nri
        vr(ir) = vsmt(i,ias)*y00
        i = i + lmmaxi
      ENDDO 
      DO ir = iro,nr
        vr(ir) = vsmt(i,ias)*y00
        i = i + lmmaxo
      ENDDO 
      !-----------------------!
      !     APW functions     !
      !-----------------------!
      DO l = 0,lmaxapw
        DO io = 1,apword(l,is)
          IF(apwve(io,l,is)) THEN 
            ! check if previous radial functions have same default energies
            DO jo = 1,io-1
              IF(apwve(jo,l,is)) THEN 
                IF(abs(apwe0(io,l,is) - apwe0(jo,l,is)) < 1.d-4) THEN 
                  apwe(io,l,ias) = apwe(jo,l,ias)
                  GOTO 10 ! next order io
                ENDIF 
              ENDIF 
            ENDDO 
            ! find the band energy starting from default
            apwe(io,l,ias) = apwe0(io,l,is)
            CALL my_findband(solsc, l, nr, rlmt(:,1,is), vr, epsband, demaxbnd, &
                          apwe(io,l,ias), fnd)
            IF(.not. fnd) nnf = nnf + 1
            !
          ELSE 
            ! set linearization energy automatically if allowed/enabled
            IF(autolinengy) apwe(io,l,ias) = efermi + dlefe
          ENDIF 
          10 CONTINUE 
        ENDDO 
      ENDDO 
      !---------------------------------!
      !     local-orbital functions     !
      !---------------------------------!
      DO ilo = 1,nlorb(is)
        DO io = 1,lorbord(ilo,is)
          IF(lorbve(io,ilo,is)) THEN 
            ! check if previous radial functions have same default energies
            DO jo = 1,io-1
              IF( lorbve(jo,ilo,is) ) THEN 
                IF(abs(lorbe0(io,ilo,is) - lorbe0(jo,ilo,is)) < 1.d-4) THEN 
                  lorbe(io,ilo,ias) = lorbe(jo,ilo,ias)
                  GOTO 20 ! next order
                ENDIF 
              ENDIF 
            ENDDO 
            l = lorbl(ilo,is)
            ! find the band energy starting from default
            lorbe(io,ilo,ias) = lorbe0(io,ilo,is)
            CALL my_findband(solsc,l,nr,rlmt(:,1,is),vr,epsband,demaxbnd, lorbe(io,ilo,ias),fnd)
            IF(.not. fnd) nnf = nnf+1
          ELSE 
            ! set linearization energy automatically if allowed/enabled
            IF(autolinengy) lorbe(io,ilo,ias) = efermi + dlefe
          ENDIF 
          20 continue ! next order
        ENDDO ! over all order
      ENDDO ! over all nlorb
      !
      done(ia) = .true.
      ! copy to equivalent atoms
      DO ja = 1,natoms(is)
        IF((.not. done(ja) ) .and. (eqatoms(ia,ja,is))) THEN 
          jas = idxas(ja,is)
          DO l = 0,lmaxapw
            DO io = 1,apword(l,is)
              apwe(io,l,jas) = apwe(io,l,ias)
            ENDDO 
          ENDDO 
          DO ilo = 1,nlorb(is)
            DO io = 1,lorbord(ilo,is)
              lorbe(io,ilo,jas) = lorbe(io,ilo,ias)
            ENDDO 
          ENDDO 
          done(ja) = .true.
        ENDIF 
      ENDDO 
  ! end loops over atoms and species
    ENDDO 
  ENDDO 
  
  IF( nnf > 0) THEN 
    WRITE(*,*)
    WRITE(*,'("Warning(linengy): could not find ",I3," linearisation energies in s.c. loop ",I5)') nnf,iscl
  ENDIF 
  
  RETURN 
END SUBROUTINE 
