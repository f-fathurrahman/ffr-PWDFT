SUBROUTINE my_genlofr()
  !USE modmain
  USE m_atoms, ONLY: natmmax, natoms, nspecies, idxas
  USE m_symmetry, ONLY: eqatoms
  USE m_apwlo, ONLY: nplorb, nlorb, lorbordmax, nlomax, lofr, lorbl, &
               lorbord, deapwlo, lorbdm, lorbe
  USE m_muffin_tins, ONLY: nrmtmax, nrmt, nrmti, rlmt, rmt, lmmaxo, lmmaxi
  USE m_constants, ONLY: y00, solsc
  USE m_density_pot_xc, ONLY: vsmt
  IMPLICIT NONE 
  ! local variables
  INTEGER :: is,ia,ja,ias,jas
  INTEGER :: nr,nri,ir,i
  INTEGER :: ilo,jlo,io,jo
  INTEGER :: nn,l,info
  REAL(8) :: e,t1
  ! automatic arrays
  LOGICAL :: done(natmmax)
  INTEGER :: ipiv(nplorb)
  REAL(8) :: vr(nrmtmax), fr(nrmtmax)
  REAL(8) :: p0(nrmtmax,lorbordmax), p1(nrmtmax)
  REAL(8) :: q0(nrmtmax), q1(nrmtmax), ep0(nrmtmax,lorbordmax)
  REAL(8) :: p0s(nrmtmax,nlomax), ep0s(nrmtmax,nlomax)
  REAL(8) :: xa(nplorb), ya(nplorb)
  REAL(8) :: a(nplorb,nplorb), b(nplorb)
  ! external functions
  REAL(8) splint, polynm
  external splint, polynm

  DO is=1,nspecies
    nr = nrmt(is)
    nri = nrmti(is)
    done(:) = .false.
    DO ia=1,natoms(is)
      IF(done(ia)) cycle
      ias=idxas(ia,is)
      ! use spherical part of potential
      i=1
      DO ir=1,nri
        vr(ir)=vsmt(i,ias)*y00
        i=i+lmmaxi
      ENDDO 
      DO ir=nri+1,nr
        vr(ir)=vsmt(i,ias)*y00
        i=i+lmmaxo
      ENDDO 
      DO ilo=1,nlorb(is)
        l=lorbl(ilo,is)
        DO jo=1,lorbord(ilo,is)
          ! linearisation energy accounting for energy derivative
          e=lorbe(jo,ilo,ias)+dble(lorbdm(jo,ilo,is))*deapwlo
          ! integrate the radial Schrodinger equation
          CALL rschrodint(solsc,l,e,nr,rlmt(:,1,is),vr,nn,p0(:,jo),p1,q0,q1)
          ep0(1:nr,jo)=e*p0(1:nr,jo)
          ! normalise radial functions
          fr(1:nr)=p0(1:nr,jo)**2
          t1=splint(nr,rlmt(:,1,is),fr)
          t1=1.d0/sqrt(abs(t1))
          CALL dscal(nr,t1,p0(:,jo),1)
          CALL dscal(nr,t1,ep0(:,jo),1)
          ! set up the matrix of radial derivatives
          DO i=1,nplorb
            ir=nr-nplorb+i
            xa(i)=rlmt(ir,1,is)
            ya(i)=p0(ir,jo)*rlmt(ir,-1,is)
          ENDDO 
          DO io=1,lorbord(ilo,is)
            a(io,jo)=polynm(io-1,nplorb,xa,ya,rmt(is))
          ENDDO 
        ENDDO 
        ! set up the target vector
        b(:)=0.d0
        b(lorbord(ilo,is))=1.d0
        CALL dgesv(lorbord(ilo,is),1,a,nplorb,ipiv,b,nplorb,info)
        IF(info.ne.0) goto 10
        ! generate linear superposition of radial functions
        p0s(:,ilo)=0.d0
        ep0s(:,ilo)=0.d0
        DO io=1,lorbord(ilo,is)
          t1=b(io)
          CALL daxpy(nr,t1,p0(:,io),1,p0s(:,ilo),1)
          CALL daxpy(nr,t1,ep0(:,io),1,ep0s(:,ilo),1)
        ENDDO 
        ! normalise radial functions
        fr(1:nr)=p0s(1:nr,ilo)**2
        t1=splint(nr,rlmt(:,1,is),fr)
        t1=1.d0/sqrt(abs(t1))
        CALL dscal(nr,t1,p0s(:,ilo),1)
        CALL dscal(nr,t1,ep0s(:,ilo),1)
        ! subtract linear combination of previous local-orbitals with same l
        DO jlo=1,ilo-1
          IF(lorbl(jlo,is).eq.l) THEN 
            fr(1:nr)=p0s(1:nr,ilo)*p0s(1:nr,jlo)
            t1=-splint(nr,rlmt(:,1,is),fr)
            CALL daxpy(nr,t1,p0s(:,jlo),1,p0s(:,ilo),1)
            CALL daxpy(nr,t1,ep0s(:,jlo),1,ep0s(:,ilo),1)
          ENDIF 
        ENDDO 
        ! normalise radial functions again
        fr(1:nr)=p0s(1:nr,ilo)**2
        t1=splint(nr,rlmt(:,1,is),fr)
        t1=abs(t1)
        IF(t1.lt.1.d-25) goto 10
        t1=1.d0/sqrt(t1)
        CALL dscal(nr,t1,p0s(:,ilo),1)
        CALL dscal(nr,t1,ep0s(:,ilo),1)
  ! divide by r and store in global array
        DO ir=1,nr
          t1=rlmt(ir,-1,is)
          lofr(ir,1,ilo,ias)=t1*p0s(ir,ilo)
          lofr(ir,2,ilo,ias)=t1*ep0s(ir,ilo)
        ENDDO 
      ENDDO 
      done(ia)=.true.
  ! copy to equivalent atoms
      DO ja=1,natoms(is)
        IF((.not.done(ja)).and.(eqatoms(ia,ja,is))) THEN 
          jas=idxas(ja,is)
          DO ilo=1,nlorb(is)
            CALL dcopy(nr,lofr(:,1,ilo,ias),1,lofr(:,1,ilo,jas),1)
            CALL dcopy(nr,lofr(:,2,ilo,ias),1,lofr(:,2,ilo,jas),1)
          ENDDO 
          done(ja)=.true.
        ENDIF 
      ENDDO 
  ! end loop over atoms and species
    ENDDO 
  ENDDO 
  RETURN 
  10 continue
  WRITE(*,*)
  WRITE(*,'("Error(genlofr): degenerate local-orbital radial functions")')
  WRITE(*,'(" for species ",I4)') is
  WRITE(*,'(" atom ",I4)') ia
  WRITE(*,'(" and local-orbital ",I4)') ilo
  WRITE(*,*)
  stop
  END SUBROUTINE 
!EOC

