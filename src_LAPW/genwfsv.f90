SUBROUTINE genwfsv(tsh,tgp,nst,idx,ngdg,igf,ngp,igpig,apwalm,evecfv,evecsv, &
 wfmt,ld,wfir)
  use modmain
  ! !INPUT/OUTPUT PARAMETERS:
  !   tsh    : .true. if wfmt should be in spherical harmonic basis (in,logical)
  !   tgp    : .true. if wfir should be in G+p-space, otherwise in real-space
  !            (in,logical)
  !   nst    : number of states to be calculated (in,integer)
  !   idx    : index to states which are to be calculated (in,integer(nst))
  !   ngdg   : G-vector grid sizes (in,integer(3))
  !   igf    : map from G-vector index to FFT array (in,integer(*))
  !   ngp    : number of G+p-vectors (in,integer(nspnfv))
  !   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv))
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  !   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
  !   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
  !   wfmt   : muffin-tin part of the wavefunctions for every state in spherical
  !            coordinates (out,complex(npcmtmax,natmtot,nspinor,nst))
  !   ld     : leading dimension of wfir (in,integer)
  !   wfir   : interstitial part of the wavefunctions for every state
  !            (out,complex(ld,nspinor,nst))
  ! !DESCRIPTION:
  !   Calculates the second-variational spinor wavefunctions in both the
  !   muffin-tin and interstitial regions for every state of a particular
  !   $k$-point. A coarse radial mesh is assumed in the muffin-tins with angular
  !   momentum cut-off of {\tt lmaxo}.
  IMPLICIT NONE 
  ! arguments
  logical, intent(in) :: tsh,tgp
  INTEGER, intent(in) :: nst,idx(nst),ngdg(3),igf(*)
  INTEGER, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
  COMPLEX(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
  COMPLEX(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
  COMPLEX(8), intent(out) :: wfmt(npcmtmax,natmtot,nspinor,nst)
  INTEGER, intent(in) :: ld
  COMPLEX(8), intent(out) :: wfir(ld,nspinor,nst)
  ! local variables
  INTEGER :: ist,ispn,jspn
  INTEGER :: is,ia,ias,i,j,k
  INTEGER :: nrc,nrci,npc
  INTEGER :: igp,ifg,nthd
  REAL(8) :: t1
  COMPLEX(8) :: zq(2),z1
  ! automatic arrays
  LOGICAL :: done(nstfv,nspnfv)
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: wfmt1(:,:,:),wfmt2(:)

  !--------------------------------!
  !     muffin-tin wavefunction    !
  !--------------------------------!
  IF(tevecsv) ALLOCATE(wfmt1(npcmtmax,nstfv,nspnfv))
  IF(.not.tsh) ALLOCATE(wfmt2(npcmtmax))
  DO ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
  ! de-phasing factor for spin-spirals
    IF(ssdph) THEN 
      ia=idxia(ias)
      t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
      zq(1)=cmplx(cos(t1),sin(t1),8)
      zq(2)=conjg(zq(1))
    ENDIF 
    done(:,:)=.false.
  ! loop only over required states
    DO j=1,nst
  ! index to state in evecsv
      k=idx(j)
      IF(tevecsv) THEN 
  ! generate spinor wavefunction from second-variational eigenvectors
        i=0
        DO ispn=1,nspinor
          jspn=jspnfv(ispn)
          wfmt(1:npc,ias,ispn,j)=0.d0
          DO ist=1,nstfv
            i=i+1
            z1=evecsv(i,k)
            IF(abs(dble(z1))+abs(aimag(z1)).gt.epsocc) THEN 
              IF(ssdph) z1=z1*zq(ispn)
              IF(.not.done(ist,jspn)) THEN 
                IF(tsh) THEN 
                  CALL wavefmt(lradstp,ias,ngp(jspn),apwalm(:,:,:,ias,jspn), &
                   evecfv(:,ist,jspn),wfmt1(:,ist,jspn))
                else
                  CALL wavefmt(lradstp,ias,ngp(jspn),apwalm(:,:,:,ias,jspn), &
                   evecfv(:,ist,jspn),wfmt2)
  ! convert to spherical coordinates
                  CALL zbsht(nrc,nrci,wfmt2,wfmt1(:,ist,jspn))
                ENDIF 
                done(ist,jspn)=.true.
              ENDIF 
  ! add to spinor wavefunction
              CALL zaxpy(npc,z1,wfmt1(:,ist,jspn),1,wfmt(:,ias,ispn,j),1)
            ENDIF 
          ENDDO 
        ENDDO 
      else
  ! spin-unpolarised wavefunction
        IF(tsh) THEN 
          CALL wavefmt(lradstp,ias,ngp,apwalm(:,:,:,ias,1),evecfv(:,k,1), &
           wfmt(:,ias,1,j))
        else
          CALL wavefmt(lradstp,ias,ngp,apwalm(:,:,:,ias,1),evecfv(:,k,1),wfmt2)
  ! convert to spherical coordinates
          CALL zbsht(nrc,nrci,wfmt2,wfmt(:,ias,1,j))
        ENDIF 
      ENDIF 
  ! end loop over second-variational states
    ENDDO 
  ! end loops over atoms
  ENDDO 
  
  IF(tevecsv) DEALLOCATE(wfmt1)
  IF(.not.tsh) DEALLOCATE(wfmt2)

  !-----------------------------------!
  !     interstitial wavefunction     !
  !-----------------------------------!
  t1=1.d0/sqrt(omega)
  DO j=1,nst
    k=idx(j)
    wfir(:,:,j)=0.d0
    IF(tevecsv) THEN 
      ! generate spinor wavefunction from second-variational eigenvectors
      i=0
      DO ispn=1,nspinor
        jspn=jspnfv(ispn)
        DO ist=1,nstfv
          i=i+1
          z1=evecsv(i,k)
          IF(abs(dble(z1))+abs(aimag(z1)).gt.epsocc) THEN 
            IF(tgp) THEN 
              ! wavefunction in G+p-space
              CALL zaxpy(ngp(jspn),z1,evecfv(:,ist,jspn),1,wfir(:,ispn,j),1)
            ELSE 
              ! wavefunction in real-space
              z1=t1*z1
              DO igp=1,ngp(jspn)
                ifg=igf(igpig(igp,jspn))
                wfir(ifg,ispn,j)=wfir(ifg,ispn,j)+z1*evecfv(igp,ist,jspn)
              ENDDO 
            ENDIF 
          ENDIF 
        ENDDO 
        ! Fourier transform wavefunction to real-space if required
        IF(.not.tgp) CALL zfftifc(3,ngdg,1,wfir(:,ispn,j))
      ENDDO 
    ELSE 
      ! spin-unpolarised wavefunction
      IF(tgp) THEN 
        CALL zcopy(ngp(1),evecfv(:,k,1),1,wfir(:,1,j),1)
      ELSE 
        DO igp=1,ngp(1)
          ifg=igf(igpig(igp,1))
          wfir(ifg,1,j)=t1*evecfv(igp,k,1)
        ENDDO 
        CALL zfftifc(3,ngdg,1,wfir(:,1,j))
      ENDIF 
    ENDIF 
  ENDDO 
  RETURN 
END SUBROUTINE 
