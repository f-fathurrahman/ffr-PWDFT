SUBROUTINE eveqnsv(ngp,igpig,vgpc,apwalm,evalfv,evecfv,evalsvp,evecsv)
  USE modmain
  use moddftu
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ngp,igpig(ngkmax)
  REAL(8), intent(in) :: vgpc(3,ngkmax)
  COMPLEX(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  REAL(8), intent(in) :: evalfv(nstfv)
  COMPLEX(8), intent(in) :: evecfv(nmatmax,nstfv)
  REAL(8), intent(out) :: evalsvp(nstsv)
  COMPLEX(8), intent(out) :: evecsv(nstsv,nstsv)
  ! local variables
  LOGICAL :: socz
  INTEGER :: nsc,nsd,ist,jst
  INTEGER :: ispn,is,ias, jspn
  INTEGER :: nrc,nrci,nrco,irco,irc
  INTEGER :: npc,npci
  INTEGER :: igp,i,j,k
  REAL(8) :: ca,t1
  REAL(8) :: ts0,ts1
  integer :: l, ld, nm, lm
  COMPLEX(8) z1
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: wfmt1(:,:),wfmt2(:),wfmt3(:),wfmt4(:,:)
  COMPLEX(8), ALLOCATABLE :: zlflm(:,:),gwfmt(:,:,:),gzfmt(:,:)
  COMPLEX(8), ALLOCATABLE :: wfir1(:),wfir2(:),wfgp(:,:),gwfgp(:,:)
  ! external functions
  COMPLEX(8) zdotc,zfmtinp
  external zdotc,zfmtinp

  ! no calculation of second-variational eigenvectors
  IF(.not.tevecsv) THEN 
    DO i=1,nstsv
      evalsvp(i)=evalfv(i)
    ENDDO 
    evecsv(:,:)=0.d0
    DO i=1,nstsv
      evecsv(i,i)=1.d0
    ENDDO 
    RETURN 
  ENDIF 
  CALL timesec(ts0)

  ! coupling constant of the external A-field (1/c)
  ca=1.d0/solsc
  
  ! number of spin combinations after application of Hamiltonian
  IF(spinpol) THEN 
    IF(ncmag.or.spinorb) THEN 
      nsc=3
    else
      nsc=2
    ENDIF 
    nsd=2
  else
    nsc=1
    nsd=1
  ENDIF 
  
  ! special case of spin-orbit coupling and collinear magnetism
  IF(spinorb.and.cmagz) THEN 
    socz=.true.
  else
    socz=.false.
  ENDIF 

  ! zero the second-variational Hamiltonian (stored in the eigenvector array)
  evecsv(:,:)=0.d0

  !-------------------------!
  !     muffin-tin part     !
  !-------------------------!
  ALLOCATE(wfmt1(npcmtmax,nstfv))
  
  IF(xcgrad.eq.4) ALLOCATE(gwfmt(npcmtmax,3,nstfv))
  
  ALLOCATE(wfmt2(npcmtmax),wfmt3(npcmtmax),wfmt4(npcmtmax,nsc))
  IF(spinorb) ALLOCATE(zlflm(lmmaxo,3))
  IF(tafield.or.(xcgrad.eq.4)) ALLOCATE(gzfmt(npcmtmax,3))

  ! begin loop over atoms
  DO ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    nrco=nrc-nrci
    irco=nrci+1
    npc=npcmt(is)
    npci=npcmti(is)
  ! compute the first-variational wavefunctions
    DO ist=1,nstfv
      CALL wavefmt(lradstp,ias,ngp,apwalm(:,:,:,ias),evecfv(:,ist),wfmt1(:,ist))
    ENDDO 
  
  ! begin loop over states
    DO jst=1,nstfv
      IF(spinpol) THEN 
  ! convert wavefunction to spherical coordinates
        CALL zbsht(nrc,nrci,wfmt1(:,jst),wfmt2)
  ! apply Kohn-Sham effective magnetic field
        wfmt3(1:npc)=bsmt(1:npc,ias,ndmag)*wfmt2(1:npc)
  ! convert to spherical harmonics and store in wfmt4
        CALL zfsht(nrc,nrci,wfmt3,wfmt4)
        wfmt4(1:npc,2)=-wfmt4(1:npc,1)
  ! non-collinear magnetic field
        IF(ncmag) THEN 
          wfmt3(1:npc)=cmplx(bsmt(1:npc,ias,1),-bsmt(1:npc,ias,2),8)*wfmt2(1:npc)
          CALL zfsht(nrc,nrci,wfmt3,wfmt4(:,3))
        ENDIF 
        IF(socz) wfmt4(1:npc,3)=0.d0
  ! apply spin-orbit coupling if required
        IF(spinorb) THEN 
  ! inner part of muffin-tin
          i=1
          DO irc=1,nrci
            CALL lopzflm(lmaxi,wfmt1(i,jst),lmmaxo,zlflm)
            t1=socfr(irc,ias)
            DO j=1,lmmaxi
              wfmt4(i,1)=wfmt4(i,1)+t1*zlflm(j,3)
              wfmt4(i,2)=wfmt4(i,2)-t1*zlflm(j,3)
              wfmt4(i,3)=wfmt4(i,3)+t1*(zlflm(j,1) &
               +cmplx(aimag(zlflm(j,2)),-dble(zlflm(j,2)),8))
              i=i+1
            ENDDO 
          ENDDO 
  ! outer part of muffin-tin
          DO irc=irco,nrc
            CALL lopzflm(lmaxo,wfmt1(i,jst),lmmaxo,zlflm)
            t1=socfr(irc,ias)
            DO j=1,lmmaxo
              wfmt4(i,1)=wfmt4(i,1)+t1*zlflm(j,3)
              wfmt4(i,2)=wfmt4(i,2)-t1*zlflm(j,3)
              wfmt4(i,3)=wfmt4(i,3)+t1*(zlflm(j,1) &
               +cmplx(aimag(zlflm(j,2)),-dble(zlflm(j,2)),8))
              i=i+1
            ENDDO 
          ENDDO 
        ENDIF 
      ELSE 
        DO k=1,nsc
          wfmt4(1:npc,k)=0.d0
        ENDDO 
      ENDIF 
  
! apply muffin-tin potential matrix if required
    if (tvmatmt) then
      do l=0,lmaxdm
        if (tvmmt(l,ias)) then
          nm=2*l+1
          lm=idxlm(l,-l)
          do k=1,nsc
            if (k.eq.1) then
              ispn=1
              jspn=1
            else if (k.eq.2) then
              ispn=2
              jspn=2
            else
              ispn=1
              jspn=2
            end if
            if (l.le.lmaxi) then
              call zgemm('N','N',nm,nrci,nm,zone,vmatmt(lm,ispn,lm,jspn,ias), &
               ld,wfmt1(lm,jst),lmmaxi,zone,wfmt4(lm,k),lmmaxi)
            end if
            i=npci+lm
            call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,ispn,lm,jspn,ias),ld, &
             wfmt1(i,jst),lmmaxo,zone,wfmt4(i,k),lmmaxo)
          end do
        end if
      end do
    end if


  ! apply vector potential if required
      IF(tafield) THEN 
        CALL gradzfmt(nrc,nrci,rlcmt(:,1,is),rlcmt(:,-1,is),wfmt1(:,jst), &
         npcmtmax,gzfmt)
        DO i=1,npc
          z1=afieldc(1)*gzfmt(i,1)+afieldc(2)*gzfmt(i,2)+afieldc(3)*gzfmt(i,3)
          z1=ca*cmplx(-aimag(z1),dble(z1),8)
          wfmt4(i,1:nsd)=wfmt4(i,1:nsd)+z1
        ENDDO 
      ENDIF 
  
  ! add to second-variational Hamiltonian matrix
  ! upper diagonal block
      DO ist=1,jst
        z1=zfmtinp(nrc,nrci,wrcmt(:,is),wfmt1(:,ist),wfmt4)
        evecsv(ist,jst)=evecsv(ist,jst)+z1
      ENDDO 
  
      ! lower diagonal block
      IF(nsc.ge.2) THEN 
        j=jst+nstfv
        DO ist=1,jst
          i=ist+nstfv
          z1=zfmtinp(nrc,nrci,wrcmt(:,is),wfmt1(:,ist),wfmt4(:,2))
          evecsv(i,j)=evecsv(i,j)+z1
        ENDDO 
      ENDIF 
  
      ! off-diagonal block
      IF(nsc.eq.3) THEN 
        DO ist=1,nstfv
          z1=zfmtinp(nrc,nrci,wrcmt(:,is),wfmt1(:,ist),wfmt4(:,3))
          evecsv(ist,j)=evecsv(ist,j)+z1
        ENDDO 
      ENDIF 
    ENDDO  ! end loop over states
  
  ! apply tau-DFT non-multiplicative potential if required
    IF(xcgrad.eq.4) THEN 
      DO ist=1,nstfv
        CALL gradzfmt(nrc,nrci,rlcmt(:,1,is),rlcmt(:,-1,is),wfmt1(:,ist), &
         npcmtmax,gwfmt(:,:,ist))
      ENDDO 
  
      DO jst=1,nstfv
        DO k=1,3
          CALL zbsht(nrc,nrci,gwfmt(:,k,jst),wfmt2)
          wfmt2(1:npc)=wsmt(1:npc,ias)*wfmt2(1:npc)
          CALL zfsht(nrc,nrci,wfmt2,gzfmt(:,k))
        ENDDO 
        DO ist=1,nstfv
          z1=0.d0
          DO k=1,3
            z1=z1+zfmtinp(nrc,nrci,wrcmt(:,is),gwfmt(:,k,ist),gzfmt(:,k))
          ENDDO 
          z1=0.5d0*z1
          evecsv(ist,jst)=evecsv(ist,jst)+z1
          IF(spinpol) THEN 
            i=ist+nstfv
            j=jst+nstfv
            evecsv(i,j)=evecsv(i,j)+z1
          ENDIF 
        ENDDO 
      ENDDO 
    ENDIF 
  
  ENDDO ! end loop over atoms
  
  DEALLOCATE(wfmt2,wfmt3,wfmt4)
  IF(spinorb) DEALLOCATE(zlflm)
  IF(tafield.or.(xcgrad.eq.4)) DEALLOCATE(gzfmt)
  
  DEALLOCATE(wfmt1)
  IF(xcgrad.eq.4) DEALLOCATE(gwfmt)

  !---------------------------!
  !     interstitial part     !
  !---------------------------!
  IF(spinpol.or.tafield.or.(xcgrad.eq.4)) THEN 
    IF(socz) nsc=2
    IF(xcgrad.eq.4) ALLOCATE(gwfgp(ngkmax,nstfv))
    ALLOCATE(wfir1(ngtot),wfir2(ngtot),wfgp(ngkmax,nsc))
  ! begin loop over states
    DO jst=1,nstfv
      wfir1(:)=0.d0
      DO igp=1,ngp
        wfir1(igfft(igpig(igp)))=evecfv(igp,jst)
      ENDDO 
  ! Fourier transform wavefunction to real-space
      CALL zfftifc(3,ngridg,1,wfir1)
  ! multiply with magnetic field and transform to G-space
      IF(spinpol) THEN 
        wfir2(:)=bsir(:,ndmag)*wfir1(:)
        CALL zfftifc(3,ngridg,-1,wfir2)
        DO igp=1,ngp
          wfgp(igp,1)=wfir2(igfft(igpig(igp)))
        ENDDO 
        wfgp(1:ngp,2)=-wfgp(1:ngp,1)
        IF(ncmag) THEN 
          wfir2(:)=cmplx(bsir(:,1),-bsir(:,2),8)*wfir1(:)
          CALL zfftifc(3,ngridg,-1,wfir2)
          DO igp=1,ngp
            wfgp(igp,3)=wfir2(igfft(igpig(igp)))
          ENDDO 
        ENDIF 
      else
        wfgp(1:ngp,1:nsd)=0.d0
      ENDIF 
  ! apply vector potential if required
      IF(tafield) THEN 
        wfir1(:)=0.d0
        DO igp=1,ngp
          t1=-ca*dot_product(afieldc(:),vgpc(:,igp))
          wfir1(igfft(igpig(igp)))=t1*evecfv(igp,jst)
        ENDDO 
        CALL zfftifc(3,ngridg,1,wfir1)
        wfir1(:)=wfir1(:)*cfunir(:)
        CALL zfftifc(3,ngridg,-1,wfir1)
        DO igp=1,ngp
          z1=wfir1(igfft(igpig(igp)))
          wfgp(igp,1:nsd)=wfgp(igp,1:nsd)+z1
        ENDDO 
      ENDIF 
  ! add to second-variational Hamiltonian matrix
  ! upper diagonal block
      DO ist=1,jst
        evecsv(ist,jst)=evecsv(ist,jst)+zdotc(ngp,evecfv(:,ist),1,wfgp,1)
      ENDDO 
  ! lower diagonal block
      IF(nsc.ge.2) THEN 
        j=jst+nstfv
        DO ist=1,jst
          i=ist+nstfv
          evecsv(i,j)=evecsv(i,j)+zdotc(ngp,evecfv(:,ist),1,wfgp(:,2),1)
        ENDDO 
      ENDIF 
  ! off-diagonal block
      IF(nsc.eq.3) THEN 
        DO ist=1,nstfv
          evecsv(ist,j)=evecsv(ist,j)+zdotc(ngp,evecfv(:,ist),1,wfgp(:,3),1)
        ENDDO 
      ENDIF 
  ! end loop over states
    ENDDO 
  
  ! apply tau-DFT non-multiplicative potential if required
    IF(xcgrad.eq.4) THEN 
      DO k=1,3
  ! determine the gradient of the wavefunctions
        DO ist=1,nstfv
          DO igp=1,ngp
            z1=evecfv(igp,ist)
            gwfgp(igp,ist)=vgpc(k,igp)*cmplx(-aimag(z1),dble(z1),8)
          ENDDO 
        ENDDO 
        DO jst=1,nstfv
          wfir1(:)=0.d0
          DO igp=1,ngp
            wfir1(igfft(igpig(igp)))=gwfgp(igp,jst)
          ENDDO 
          CALL zfftifc(3,ngridg,1,wfir1)
          wfir1(:)=wsir(:)*wfir1(:)
          CALL zfftifc(3,ngridg,-1,wfir1)
          DO igp=1,ngp
            wfgp(igp,1)=wfir1(igfft(igpig(igp)))
          ENDDO 
          DO ist=1,nstfv
            z1=0.5d0*zdotc(ngp,gwfgp(:,ist),1,wfgp,1)
            evecsv(ist,jst)=evecsv(ist,jst)+z1
            IF(spinpol) THEN 
              i=ist+nstfv
              j=jst+nstfv
              evecsv(i,j)=evecsv(i,j)+z1
            ENDIF 
          ENDDO 
        ENDDO 
      ENDDO 
    ENDIF 
    DEALLOCATE(wfir1,wfir2,wfgp)
    IF(xcgrad.eq.4) DEALLOCATE(gwfgp)
  ENDIF 

  ! add the diagonal first-variational part
  i=0
  DO ispn=1,nspinor
    DO ist=1,nstfv
      i=i+1
      evecsv(i,i)=evecsv(i,i)+evalfv(ist)
    ENDDO 
  ENDDO 

  IF(spcpl.or.(.not.spinpol)) THEN 
  ! spins are coupled; or spin-unpolarised: full diagonalisation
    CALL eveqnz(nstsv,nstsv,evecsv,evalsvp)
  ELSE 
  ! spins not coupled: block diagonalise H
    CALL eveqnz(nstfv,nstsv,evecsv,evalsvp)
    i=nstfv+1
    CALL eveqnz(nstfv,nstsv,evecsv(i,i),evalsvp(i))
    DO i=1,nstfv
      DO j=1,nstfv
        evecsv(i,j+nstfv)=0.d0
        evecsv(i+nstfv,j)=0.d0
      ENDDO 
    ENDDO 
  ENDIF 

  CALL timesec(ts1)
  timesv=timesv+ts1-ts0
  
  RETURN 
END SUBROUTINE 

