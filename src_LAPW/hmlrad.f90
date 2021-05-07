SUBROUTINE hmlrad()
  USE m_atoms, ONLY: idxis, natmtot
  USE m_muffin_tins, ONLY: nrmt, nrmti, nrmtmax, npmti, lmmaxi, lmmaxo, lmaxo, lmaxi, &
                   lmaxapw, wrmt
  USE m_apwlo, ONLY: nlorb, lofr, lorbl, apword, apwfr
  USE m_hamiltonian, ONLY: hlolo, hloa, haa
  USE m_density_pot_xc, ONLY: vsmt
  USE m_constants, ONLY: y00
  IMPLICIT NONE 
  ! local variables
  INTEGER is,ias
  INTEGER nr,nri,iro
  INTEGER ir,npi,i
  INTEGER l1,l2,l3,m2,lm2
  INTEGER io,jo,ilo,jlo
  REAL(8) t1
  ! allocatable arrays 
  REAL(8), ALLOCATABLE :: fr(:)
  ! begin loops over atoms and species
  ALLOCATE( fr(nrmtmax) )
  
  DO ias = 1,natmtot
    
    is = idxis(ias)
    nr = nrmt(is)
    nri = nrmti(is)
    iro = nri + 1
    npi = npmti(is)

    !---------------------------!
    !     APW-APW integrals     !
    !---------------------------!
    DO l1=0,lmaxapw
      DO io=1,apword(l1,is)
        DO l3=0,lmaxapw
          DO jo=1,apword(l3,is)
            IF( l1 == l3 ) THEN 
              fr(1:nr)=apwfr(1:nr,1,io,l1,ias)*apwfr(1:nr,2,jo,l3,ias)
              t1=dot_product(wrmt(1:nr,is),fr(1:nr))
              haa(1,jo,l3,io,l1,ias)=t1/y00
            ELSE 
              haa(1,jo,l3,io,l1,ias)=0.d0
            ENDIF 
            IF( l1 >= l3 ) THEN 
              lm2=1
              DO l2=1,lmaxi
                DO m2=-l2,l2
                  lm2=lm2+1
                  i=lm2
                  DO ir=1,nri
                    t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)
                    fr(ir)=t1*vsmt(i,ias)
                    i=i+lmmaxi
                  ENDDO 
                  DO ir=iro,nr
                    t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)
                    fr(ir)=t1*vsmt(i,ias)
                    i=i+lmmaxo
                  ENDDO 
                  t1=dot_product(wrmt(1:nr,is),fr(1:nr))
                  haa(lm2,jo,l3,io,l1,ias)=t1
                  haa(lm2,io,l1,jo,l3,ias)=t1
                ENDDO 
              ENDDO 
              DO l2=lmaxi+1,lmaxo
                DO m2=-l2,l2
                  lm2=lm2+1
                  i=npi+lm2
                  DO ir=iro,nr
                    t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)
                    fr(ir)=t1*vsmt(i,ias)
                    i=i+lmmaxo
                  ENDDO 
                  t1=dot_product(wrmt(iro:nr,is),fr(iro:nr))
                  haa(lm2,jo,l3,io,l1,ias)=t1
                  haa(lm2,io,l1,jo,l3,ias)=t1
                ENDDO 
              ENDDO 
            ENDIF 
          ENDDO 
        ENDDO 
      ENDDO 
    ENDDO 
  
    !-------------------------------------!
    !     local-orbital-APW integrals     !
    !-------------------------------------!
    DO ilo = 1,nlorb(is)
      l1 = lorbl(ilo,is)
      DO l3 = 0,lmaxapw
        DO io = 1,apword(l3,is)
          IF( l1 == l3 ) THEN 
            fr(1:nr) = lofr(1:nr,1,ilo,ias)*apwfr(1:nr,2,io,l3,ias)
            t1 = dot_product(wrmt(1:nr,is),fr(1:nr))
            hloa(1,io,l3,ilo,ias) = t1/y00
          ELSE 
            hloa(1,io,l3,ilo,ias) = 0.d0
          ENDIF 
          lm2 = 1
          DO l2=1,lmaxi
            DO m2=-l2,l2
              lm2 = lm2+1
              i = lm2
              DO ir=1,nri
                t1 = lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)
                fr(ir) = t1*vsmt(i,ias)
                i = i+lmmaxi
              ENDDO 
              DO ir=nri+1,nr
                t1 = lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)
                fr(ir) = t1*vsmt(i,ias)
                i = i + lmmaxo
              ENDDO 
              t1 = dot_product(wrmt(1:nr,is),fr(1:nr))
              hloa(lm2,io,l3,ilo,ias) = t1
            ENDDO 
          ENDDO 
          DO l2=lmaxi+1,lmaxo
            DO m2=-l2,l2
              lm2=lm2+1
              i=npi+lm2
              DO ir=iro,nr
                t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)
                fr(ir)=t1*vsmt(i,ias)
                i=i+lmmaxo
              ENDDO 
              t1=dot_product(wrmt(iro:nr,is),fr(iro:nr))
              hloa(lm2,io,l3,ilo,ias)=t1
            ENDDO 
          enddo
        ENDDO 
      ENDDO 
    ENDDO 
  
    !-----------------------------------------------!
    !     local-orbital-local-orbital integrals     !
    !-----------------------------------------------!
    DO ilo=1,nlorb(is)
      l1=lorbl(ilo,is)
      DO jlo=1,nlorb(is)
        l3=lorbl(jlo,is)
        IF( l1 == l3 ) THEN 
          fr(1:nr)=lofr(1:nr,1,ilo,ias)*lofr(1:nr,2,jlo,ias)
          t1=dot_product(wrmt(1:nr,is),fr(1:nr))
          hlolo(1,jlo,ilo,ias)=t1/y00
        ELSE 
          hlolo(1,jlo,ilo,ias)=0.d0
        ENDIF 
        lm2=1
        DO l2=1,lmaxi
          DO m2=-l2,l2
            lm2=lm2+1
            i=lm2
            DO ir=1,nri
              t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)
              fr(ir)=t1*vsmt(i,ias)
              i=i+lmmaxi
            ENDDO 
            DO ir=iro,nr
              t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)
              fr(ir)=t1*vsmt(i,ias)
              i=i+lmmaxo
            ENDDO 
            t1=dot_product(wrmt(1:nr,is),fr(1:nr))
            hlolo(lm2,jlo,ilo,ias)=t1
          ENDDO 
        ENDDO 
        DO l2=lmaxi+1,lmaxo
          DO m2=-l2,l2
            lm2 = lm2+1
            i = npi + lm2
            DO ir=iro,nr
              t1 = lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)
              fr(ir) = t1*vsmt(i,ias)
              i = i+lmmaxo
            ENDDO 
            t1 = dot_product(wrmt(iro:nr,is),fr(iro:nr))
            hlolo(lm2,jlo,ilo,ias) = t1
          ENDDO 
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO     ! end loops over atoms and species
  DEALLOCATE(fr)
  RETURN 
END SUBROUTINE 