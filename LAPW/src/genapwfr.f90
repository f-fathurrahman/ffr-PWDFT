SUBROUTINE genapwfr()
  USE m_atoms, ONLY: natmmax, natoms, idxas, nspecies
  USE m_symmetry, ONLY: eqatoms
  USE m_muffin_tins, ONLY: rlmt, nrmt, nrmti, nrmtmax, lmmaxi, lmmaxo, lmaxapw
  USE m_apwlo, ONLY: apword, apwfr, apwordmax, apwdfr, deapwlo, apwdm, apwe
  USE m_density_pot_xc, ONLY: vsmt
  USE m_constants, ONLY: solsc, y00
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ia,ja,ias,jas
  INTEGER :: nr,nri,ir,i
  INTEGER :: nn,l,io,jo
  REAL(8) :: e,t1
  ! automatic arrays
  LOGICAL :: done(natmmax)
  REAL(8) :: vr(nrmtmax),fr(nrmtmax)
  REAL(8) :: p0(nrmtmax,apwordmax),p1(nrmtmax),p1s(apwordmax)
  REAL(8) :: q0(nrmtmax),q1(nrmtmax),ep0(nrmtmax,apwordmax)
  ! external functions
  REAL(8) splint
  EXTERNAL splint

  DO is=1,nspecies
    nr = nrmt(is)
    nri = nrmti(is)
    done(:) = .false.
    DO ia=1,natoms(is)
      IF( done(ia) ) CYCLE 
      ias = idxas(ia,is)
      ! use spherical part of potential
      i = 1
      DO ir = 1,nri
        vr(ir) = vsmt(i,ias)*y00
        i = i + lmmaxi
      ENDDO 
      !
      DO ir=nri+1,nr
        vr(ir) = vsmt(i,ias)*y00
        i = i + lmmaxo
      ENDDO 
      !
      DO l=0,lmaxapw
        DO io=1,apword(l,is)
          ! linearisation energy accounting for energy derivative
          e = apwe(io,l,ias) + dble(apwdm(io,l,is))*deapwlo
          ! integrate the radial Schrodinger equation
          CALL rschrodint(solsc,l,e,nr,rlmt(:,1,is),vr,nn,p0(:,io),p1,q0,q1)
          ! multiply by the linearisation energy
          ep0(1:nr,io) = e*p0(1:nr,io)
          ! normalise radial functions
          fr(1:nr) = p0(1:nr,io)**2
          t1 = splint(nr,rlmt(:,1,is),fr)
          t1 = 1.d0/sqrt(abs(t1))
          CALL dscal(nr,t1,p0(:,io),1)
          p1s(io) = t1*p1(nr)
          CALL dscal(nr,t1,ep0(:,io),1)
          ! subtract linear combination of previous vectors
          DO jo=1,io-1
            fr(1:nr) = p0(1:nr,io)*p0(1:nr,jo)
            t1 = -splint(nr,rlmt(:,1,is),fr)
            CALL daxpy(nr,t1,p0(:,jo),1,p0(:,io),1)
            p1s(io) = p1s(io) + t1*p1s(jo)
            CALL daxpy(nr, t1, ep0(:,jo), 1, ep0(:,io), 1)
          ENDDO 
          ! normalise radial functions again
          fr(1:nr) = p0(1:nr,io)**2
          t1 = splint(nr,rlmt(:,1,is),fr)
          t1 = abs(t1)
          IF( t1 < 1.d-25) THEN 
            WRITE(*,*)
            WRITE(*,'("Error(genapwfr): degenerate APW radial functions")')
            WRITE(*,'(" for species ",I4)') is
            WRITE(*,'(" atom ",I4)') ia
            WRITE(*,'(" angular momentum ",I4)') l
            WRITE(*,'(" and order ",I4)') io
            WRITE(*,*)
            STOP 
          ENDIF 
          t1 = 1.d0/sqrt(t1)
          CALL dscal(nr,t1,p0(:,io),1)
          p1s(io) = t1*p1s(io)
          CALL dscal(nr,t1,ep0(:,io),1)
          ! divide by r and store in global array
          DO ir=1,nr
            t1=rlmt(ir,-1,is)
            apwfr(ir,1,io,l,ias)=t1*p0(ir,io)
            apwfr(ir,2,io,l,ias)=t1*ep0(ir,io)
          ENDDO 
          ! derivative at the muffin-tin surface
          apwdfr(io,l,ias)=(p1s(io)-p0(nr,io)*t1)*t1
        ENDDO 
      ENDDO 
      done(ia) = .true.
      ! copy to equivalent atoms
      DO ja=1,natoms(is)
        IF( (.not. done(ja) ) .and. ( eqatoms(ia,ja,is) ) ) THEN 
          jas=idxas(ja,is)
          DO l=0,lmaxapw
            DO io=1,apword(l,is)
              CALL dcopy(nr,apwfr(:,1,io,l,ias),1,apwfr(:,1,io,l,jas),1)
              CALL dcopy(nr,apwfr(:,2,io,l,ias),1,apwfr(:,2,io,l,jas),1)
              apwdfr(io,l,jas) = apwdfr(io,l,ias)
            ENDDO 
          ENDDO 
          done(ja)=.true.
        ENDIF 
      ENDDO 
      ! end loop over atoms and species
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 


