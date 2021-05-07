SUBROUTINE eveqnss(ngp,igpig,apwalm,evalfv,evecfv,evalsvp,evecsv)
  USE modmain
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
  COMPLEX(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
  REAL(8), intent(in) :: evalfv(nstfv,nspnfv)
  COMPLEX(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
  REAL(8), intent(out) :: evalsvp(nstsv)
  COMPLEX(8), intent(out) :: evecsv(nstsv,nstsv)
  ! local variables
  INTEGER :: ist,jst,ispn,jspn
  INTEGER :: is,ia,ias,i,j,k
  INTEGER :: nrc,nrci,nrco
  INTEGER :: l,lm,nm,npc,npci
  INTEGER :: igp,nthd
  REAL(8) :: t1
  REAL(8) :: ts0,ts1
  COMPLEX(8) zq,z1
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: wfmt1(:,:,:),wfmt2(:,:),wfmt3(:),wfmt4(:,:)
  COMPLEX(8), ALLOCATABLE :: wfir1(:,:),wfir2(:),wfgp(:,:)
  ! external functions
  COMPLEX(8) zdotc,zfmtinp
  external zdotc,zfmtinp
  
  IF(.not.spinpol) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(eveqnss): spin-unpolarised calculation")')
    WRITE(*,*)
    stop
  ENDIF 

  CALL timesec(ts0)

  ! zero the second-variational Hamiltonian (stored in the eigenvector array)
  evecsv(:,:)=0.d0
  
  !-------------------------!
  !     muffin-tin part     !
  !-------------------------!
  ALLOCATE(wfmt1(npcmtmax,nstfv,nspnfv))
  ALLOCATE(wfmt2(npcmtmax,nspnfv),wfmt3(npcmtmax),wfmt4(npcmtmax,3))
  DO ias=1,natmtot
    is=idxis(ias)
    ia=idxia(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    nrco=nrc-nrci
    npc=npcmt(is)
    npci=npcmti(is)
    ! de-phasing factor (FC, FB & LN)
    t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
    zq=cmplx(cos(t1),sin(t1),8)
    ! compute the first-variational wavefunctions
    DO ispn=1,nspnfv
      IF(ispn.eq.2) zq=conjg(zq)
      DO ist=1,nstfv
        CALL wavefmt(lradstp,ias,ngp(ispn),apwalm(:,:,:,ias,ispn), &
         evecfv(:,ist,ispn),wfmt1(:,ist,ispn))
        ! de-phase if required
        IF(ssdph) wfmt1(1:npc,ist,ispn)=zq*wfmt1(1:npc,ist,ispn)
      ENDDO 
    ENDDO 
    DO jst=1,nstfv
      ! convert wavefunction to spherical coordinates
      DO ispn=1,nspnfv
        CALL zbsht(nrc,nrci,wfmt1(:,jst,ispn),wfmt2(:,ispn))
      ENDDO 
      ! apply effective magnetic field and convert to spherical harmonics
      wfmt3(1:npc)=bsmt(1:npc,ias,3)*wfmt2(1:npc,1)
      CALL zfsht(nrc,nrci,wfmt3,wfmt4)
      wfmt3(1:npc)=-bsmt(1:npc,ias,3)*wfmt2(1:npc,2)
      CALL zfsht(nrc,nrci,wfmt3,wfmt4(:,2))
      wfmt3(1:npc)=cmplx(bsmt(1:npc,ias,1),-bsmt(1:npc,ias,2),8)*wfmt2(1:npc,2)
      CALL zfsht(nrc,nrci,wfmt3,wfmt4(:,3))
      
      ! add to second-variational Hamiltonian matrix
      ! upper diagonal block
      DO ist=1,jst
        z1=zfmtinp(nrc,nrci,wrcmt(:,is),wfmt1(:,ist,1),wfmt4)
        evecsv(ist,jst)=evecsv(ist,jst)+z1
      ENDDO 
      ! lower diagonal block
      j=jst+nstfv
      DO ist=1,jst
        i=ist+nstfv
        z1=zfmtinp(nrc,nrci,wrcmt(:,is),wfmt1(:,ist,2),wfmt4(:,2))
        evecsv(i,j)=evecsv(i,j)+z1
      ENDDO 
      ! off-diagonal block
      DO ist=1,nstfv
        z1=zfmtinp(nrc,nrci,wrcmt(:,is),wfmt1(:,ist,1),wfmt4(:,3))
        evecsv(ist,j)=evecsv(ist,j)+z1
      ENDDO 
    ENDDO 
  ! end loop over atoms
  ENDDO 
  DEALLOCATE(wfmt2,wfmt3,wfmt4)
  DEALLOCATE(wfmt1)
  
  !---------------------------!
  !     interstitial part     !
  !---------------------------!
  ALLOCATE(wfir1(ngtot,nspnfv),wfir2(ngtot),wfgp(ngkmax,3))
  ! begin loop over states
  DO jst=1,nstfv
    DO ispn=1,nspnfv
      wfir1(:,ispn)=0.d0
      DO igp=1,ngp(ispn)
        wfir1(igfft(igpig(igp,ispn)),ispn)=evecfv(igp,jst,ispn)
      ENDDO 
  ! Fourier transform wavefunction to real-space
      CALL zfftifc(3,ngridg,1,wfir1(:,ispn))
    ENDDO 
  ! multiply with magnetic field and transform to G-space
    wfir2(:)=bsir(:,3)*wfir1(:,1)
    CALL zfftifc(3,ngridg,-1,wfir2)
    DO igp=1,ngp(1)
      wfgp(igp,1)=wfir2(igfft(igpig(igp,1)))
    ENDDO 
    wfir2(:)=-bsir(:,3)*wfir1(:,2)
    CALL zfftifc(3,ngridg,-1,wfir2)
    DO igp=1,ngp(2)
      wfgp(igp,2)=wfir2(igfft(igpig(igp,2)))
    ENDDO 
    wfir2(:)=cmplx(bsir(:,1),-bsir(:,2),8)*wfir1(:,2)
    CALL zfftifc(3,ngridg,-1,wfir2)
    DO igp=1,ngp(1)
      wfgp(igp,3)=wfir2(igfft(igpig(igp,1)))
    ENDDO 
  ! add to second-variational Hamiltonian matrix
  ! upper diagonal block
    DO ist=1,jst
      evecsv(ist,jst)=evecsv(ist,jst)+zdotc(ngp(1),evecfv(:,ist,1),1,wfgp(:,1),1)
    ENDDO 
  ! lower diagonal block
    j=jst+nstfv
    DO ist=1,jst
      i=ist+nstfv
      evecsv(i,j)=evecsv(i,j)+zdotc(ngp(2),evecfv(:,ist,2),1,wfgp(:,2),1)
    ENDDO 
  ! off-diagonal block
    DO ist=1,nstfv
      evecsv(ist,j)=evecsv(ist,j)+zdotc(ngp(1),evecfv(:,ist,1),1,wfgp(:,3),1)
    ENDDO 
  ENDDO 
  DEALLOCATE(wfir1,wfir2,wfgp)

  ! add the diagonal first-variational part
  i=0
  DO ispn=1,nspinor
    DO ist=1,nstfv
      i=i+1
      evecsv(i,i)=evecsv(i,i)+evalfv(ist,ispn)
    ENDDO 
  ENDDO 
  ! diagonalise the second-variational Hamiltonian
  CALL eveqnz(nstsv,nstsv,evecsv,evalsvp)
  CALL timesec(ts1)
  timesv=timesv+ts1-ts0
  RETURN 
END SUBROUTINE 

