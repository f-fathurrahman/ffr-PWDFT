SUBROUTINE gencore()
  USE m_constants, ONLY: y00, solsc
  USE m_atoms, ONLY: natmmax, nrspmax, natoms, nstspmax, nstsp, spcore, idxas, &
                 vrsp, nrspmax, nrsp, nspecies, rsp, nsp, lsp, ksp
  USE m_symmetry, ONLY: eqatoms
  USE m_muffin_tins, ONLY: nrmt, nrmtmax, nrmti, rlmt
  USE m_spin, ONLY: ndmag, ncmag
  USE m_density_pot_xc, ONLY: vsmt, bxcmt
  USE m_core_states, ONLY: evalcr, rhocr, rwfcr, occcr, nspncr, spincore
  IMPLICIT NONE 
  ! local variables
  INTEGER :: ist,ispn,idm
  INTEGER :: is,ia,ja,ias,jas
  INTEGER :: nr,nri,nrs,ir
  REAL(8) :: t1
  ! automatic arrays
  logical done(natmmax)
  REAL(8) vr(nrspmax),eval(nstspmax)
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: br(:),fr(:,:)

  IF (spincore) ALLOCATE(br(nrmtmax),fr(nrmtmax,ndmag))
  
  ! loop over species and atoms
  DO is=1,nspecies
    nr=nrmt(is)
    nri=nrmti(is)
    nrs=nrsp(is)
    done(:)=.false.
    DO ia=1,natoms(is)
      IF( done(ia) ) CYCLE 
      ias = idxas(ia,is)
      ! Kohn-Sham magnetic field for spin-polarised core
      IF( spincore ) THEN 
        DO idm=1,ndmag
          CALL rfmtlm(1,nr,nri,bxcmt(:,ias,idm),fr(:,idm))
        ENDDO 
        IF( ncmag ) THEN 
          DO ir=1,nr
            br(ir)=sqrt(fr(ir,1)**2+fr(ir,2)**2+fr(ir,3)**2)*y00
          ENDDO 
        ELSE
          DO ir=1,nr
            br(ir)=abs(fr(ir,1))*y00
          ENDDO 
        ENDIF 
      ENDIF 
      !
      ! loop over spin channels
      DO ispn=1,nspncr
        ! use the spherical part of the crystal Kohn-Sham potential
        CALL rfmtlm(1, nr, nri, vsmt(:,ias), vr)
        vr(1:nr)=vr(1:nr)*y00
        ! spin-up and -down potentials for polarised core
        IF(spincore) THEN 
          IF(ispn.eq.1) THEN 
            vr(1:nr)=vr(1:nr)+br(1:nr)
          else
            vr(1:nr)=vr(1:nr)-br(1:nr)
          ENDIF 
        ENDIF 
        ! append the Kohn-Sham potential from the atomic calculation for r > R_MT
        t1=vr(nr)-vrsp(nr,is)
        DO ir=nr+1,nrs
          vr(ir)=vrsp(ir,is)+t1
        ENDDO 
        rhocr(:,ias,ispn)=0.d0
        DO ist=1,nstsp(is)
          IF( spcore(ist,is) ) THEN 
            ! solve the Dirac equation
            eval(ist)=evalcr(ist,ias)
            CALL rdirac(solsc,nsp(ist,is),lsp(ist,is),ksp(ist,is),nrs,rsp(:,is), &
             vr,eval(ist),rwfcr(:,1,ist,ias),rwfcr(:,2,ist,ias))
            IF( spincore ) THEN 
              ! use the spin-averaged eigenvalue for the polarised core
              IF(ispn.eq.1) THEN 
                evalcr(ist,ias)=eval(ist)
              else
                evalcr(ist,ias)=0.5d0*(evalcr(ist,ias)+eval(ist))
              ENDIF 
              t1=0.5d0*occcr(ist,ias)
            else
              evalcr(ist,ias)=eval(ist)
              t1=occcr(ist,ias)
            ENDIF 
            ! add to the core density
            DO ir=1,nr
              rhocr(ir,ias,ispn)=rhocr(ir,ias,ispn) &
               +t1*(rwfcr(ir,1,ist,ias)**2+rwfcr(ir,2,ist,ias)**2)
            ENDDO 
          ENDIF 
        ENDDO 
        DO ir=1,nr
          rhocr(ir,ias,ispn)=rhocr(ir,ias,ispn)*rlmt(ir,-2,is)*y00
        ENDDO 
      ! end loop over spin channels
      ENDDO 
      done(ia)=.true.
      ! copy to equivalent atoms
      DO ja=1,natoms(is)
        IF((.not.done(ja)).and.(eqatoms(ia,ja,is))) THEN 
          jas=idxas(ja,is)
          DO ist=1,nstsp(is)
            IF(spcore(ist,is)) THEN 
              evalcr(ist,jas)=evalcr(ist,ias)
              rwfcr(1:nrs,:,ist,jas)=rwfcr(1:nrs,:,ist,ias)
            ENDIF 
          ENDDO 
          rhocr(1:nr,jas,:)=rhocr(1:nr,ias,:)
          done(ja)=.true.
        ENDIF 
      ENDDO 
    ENDDO ! end loop over species and atoms
  ENDDO 
  IF(spincore) DEALLOCATE(br,fr)
  RETURN 
END SUBROUTINE 