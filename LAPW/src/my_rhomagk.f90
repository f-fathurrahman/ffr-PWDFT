
! input information: ik, pwgrid, apwalm, evecfv, and evecsv
!---------------------------------------------------------------------
SUBROUTINE my_rhomagk(ngp,igpig,wppt,occsvp,apwalm,evecfv,evecsv)
!---------------------------------------------------------------------
  USE m_gvectors, ONLY: igfft, ngtot, ngridg
  USE m_gkvectors, ONLY: ngkmax
  USE m_spin, ONLY: jspnfv, ssdph, spinpol, nspinor, ncmag, nspnfv, vqcss
  USE m_atoms, ONLY: idxia, idxis, natmtot, atposc
  USE m_muffin_tins, ONLY: npcmt, nrcmt, nrcmti, npcmt, npcmtmax, lradstp, &
     & lmmaxapw
  USE m_states, ONLY: tevecsv, epsocc, nstsv, nstfv
  USE m_apwlo, ONLY: apwordmax
  USE m_hamiltonian, ONLY: nmatmax
  USE m_density_pot_xc, ONLY: rhoir, rhomt, magir, magmt
  USE m_lattice, ONLY: omega
  !
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ngp(nspnfv), igpig(ngkmax,nspnfv)
  REAL(8), intent(in) :: wppt, occsvp(nstsv)
  COMPLEX(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
  COMPLEX(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
  ! local variables
  INTEGER :: ispn,jspn,ist
  INTEGER :: is,ia,ias,i,j
  INTEGER :: nrc,nrci,npc
  INTEGER :: igp,ifg
  REAL(8) :: wo,t1
  COMPLEX(8) :: zq(2),z1
  ! automatic arrays
  LOGICAL :: done(nstfv,nspnfv)
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: wfmt1(:,:,:),wfmt2(:),wfmt3(:,:),wfir(:,:)

  !----------------------------------------------!
  !     muffin-tin density and magnetisation     !
  !----------------------------------------------!
  IF(tevecsv) ALLOCATE( wfmt1(npcmtmax,nstfv,nspnfv) )
  ALLOCATE( wfmt2(npcmtmax), wfmt3(npcmtmax,nspinor) )
  !
  ! loop over all atoms
  DO ias = 1,natmtot
    is = idxis(ias)
    nrc = nrcmt(is)
    nrci = nrcmti(is)
    npc = npcmt(is)
    ! de-phasing factor for spin-spirals
    IF(ssdph) THEN 
      ia = idxia(ias)
      t1 = -0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
      zq(1) = cmplx(cos(t1),sin(t1),8)
      zq(2) = conjg(zq(1))
    ENDIF 
    
    done(:,:) = .false.
    ! loop over all second-variational state
    DO j = 1,nstsv
      wo = occsvp(j) ! occupation number
      ! skip this state if it is empty or nearly empty
      IF(abs(wo) < epsocc) cycle
      !
      wo = wo*wppt
      !
      ! if using 2nd variational scheme
      IF( tevecsv ) THEN 
        ! generate spinor wavefunction from second-variational eigenvectors
        i = 0
        DO ispn = 1,nspinor
          jspn = jspnfv(ispn)
          wfmt3(1:npc,ispn) = 0.d0
          DO ist = 1,nstfv
            i = i + 1
            z1 = evecsv(i,j)
            IF(abs(dble(z1))+abs(aimag(z1)).gt.epsocc) THEN 
              IF(ssdph) z1=z1*zq(ispn)
              IF(.not.done(ist,jspn)) THEN 
                CALL wavefmt(lradstp, ias, ngp(jspn), apwalm(:,:,:,ias,jspn), evecfv(:,ist,jspn), wfmt2)
                ! convert to spherical coordinates
                CALL zbsht(nrc,nrci,wfmt2,wfmt1(:,ist,jspn))
                done(ist,jspn)=.true.
              ENDIF 
              ! add to spinor wavefunction
              CALL zaxpy(npc,z1,wfmt1(:,ist,jspn),1,wfmt3(:,ispn),1)
            ENDIF 
          ENDDO 
        ENDDO 
        !
      ELSE 
        ! not using 2nd variational scheme
        !
        ! spin-unpolarised wavefunction
        CALL wavefmt(lradstp, ias, ngp, apwalm(:,:,:,ias,1), evecfv(:,j,1), wfmt2)
        ! convert to spherical coordinates
        CALL zbsht(nrc, nrci, wfmt2, wfmt3)
      ENDIF 
      !
      ! add to density and magnetisation
      IF(spinpol) THEN 
        ! spin-polarised
        IF(ncmag) THEN 
          ! non-collinear
          CALL rmk1(npc,wo,wfmt3,wfmt3(:,2),rhomt(:,ias),magmt(:,ias,1), &
           magmt(:,ias,2),magmt(:,ias,3))
        ELSE 
          ! collinear
          CALL rmk2(npc,wo,wfmt3,wfmt3(:,2),rhomt(:,ias),magmt(:,ias,1))
        ENDIF 
      ELSE
        ! spin-unpolarised
        CALL rmk3(npc, wo, wfmt3, rhomt(:,ias))
      ENDIF 
  
    ENDDO ! over states
  
  ENDDO ! end loop over atoms

  IF(tevecsv) DEALLOCATE(wfmt1)
  DEALLOCATE(wfmt2,wfmt3)


  !------------------------------------------------!
  !     interstitial density and magnetisation     !
  !------------------------------------------------!
  ALLOCATE(wfir(ngtot,nspinor))
  !
  ! loop over all states
  !
  DO j = 1,nstsv
    wo = occsvp(j)
    !
    IF(abs(wo) < epsocc) cycle
    !
    wo = wo*wppt/omega
    wfir(:,:) = 0.d0
    !
    IF(tevecsv) THEN 
      ! generate spinor wavefunction from second-variational eigenvectors
      i=0
      DO ispn=1,nspinor
        jspn=jspnfv(ispn)
        DO ist=1,nstfv
          i=i+1
          z1=evecsv(i,j)
          IF(abs(dble(z1))+abs(aimag(z1)).gt.epsocc) THEN 
            DO igp=1,ngp(jspn)
              ifg=igfft(igpig(igp,jspn))
              wfir(ifg,ispn)=wfir(ifg,ispn)+z1*evecfv(igp,ist,jspn)
            ENDDO 
          ENDIF 
        ENDDO 
      ENDDO 
    ELSE 
      ! spin-unpolarised wavefunction
      DO igp = 1,ngp(1)
        ifg = igfft(igpig(igp,1))
        wfir(ifg,1) = evecfv(igp,j,1)
      ENDDO 
    ENDIF 
    
    ! Fourier transform wavefunction to real-space
    DO ispn = 1,nspinor
      CALL zfftifc(3,ngridg,1,wfir(:,ispn))
    ENDDO 
    ! add to density and magnetisation
    IF( spinpol ) THEN 
      ! spin-polarised
      IF(ncmag) THEN 
        ! non-collinear
        CALL rmk1(ngtot,wo,wfir,wfir(:,2),rhoir,magir,magir(:,2),magir(:,3))
      ELSE 
        ! collinear
        CALL rmk2(ngtot,wo,wfir,wfir(:,2),rhoir,magir)
      ENDIF 
    ELSE 
      ! spin-unpolarised
      CALL rmk3(ngtot, wo, wfir, rhoir)
    ENDIF 
  ENDDO ! nstsv

  DEALLOCATE(wfir)
  RETURN 

CONTAINS 

!----------------------------------------------------
PURE SUBROUTINE rmk1(n,wo,wf1,wf2,rho,mag1,mag2,mag3)
!---------------------------------------------------
IMPLICIT NONE 
! arguments
INTEGER, intent(in) :: n
REAL(8), intent(in) :: wo
COMPLEX(8), intent(in) :: wf1(n),wf2(n)
REAL(8), intent(inout) :: rho(n),mag1(n),mag2(n),mag3(n)
! local variables
INTEGER i
REAL(8) wo2,t1,t2
COMPLEX(8) z1,z2
wo2=2.d0*wo
DO i=1,n
  z1=wf1(i); z2=wf2(i)
  t1=dble(z1)**2+aimag(z1)**2
  t2=dble(z2)**2+aimag(z2)**2
  z1=conjg(z1)*z2
  mag1(i)=mag1(i)+wo2*dble(z1)
  mag2(i)=mag2(i)+wo2*aimag(z1)
  mag3(i)=mag3(i)+wo*(t1-t2)
  rho(i)=rho(i)+wo*(t1+t2)
ENDDO 
RETURN 
END SUBROUTINE 

PURE SUBROUTINE rmk2(n,wo,wf1,wf2,rho,mag)
IMPLICIT NONE 
! arguments
INTEGER, intent(in) :: n
REAL(8), intent(in) :: wo
COMPLEX(8), intent(in) :: wf1(n),wf2(n)
REAL(8), intent(inout) :: rho(n),mag(n)
! local variables
INTEGER i
REAL(8) t1,t2
DO i=1,n
  t1=dble(wf1(i))**2+aimag(wf1(i))**2
  t2=dble(wf2(i))**2+aimag(wf2(i))**2
  mag(i)=mag(i)+wo*(t1-t2)
  rho(i)=rho(i)+wo*(t1+t2)
ENDDO 
RETURN 
END SUBROUTINE 


! used for non-spin-polarized case
PURE SUBROUTINE rmk3(n,wo,wf,rho)
IMPLICIT NONE 
! arguments
INTEGER, intent(in) :: n
REAL(8), intent(in) :: wo
COMPLEX(8), intent(in) :: wf(n)
REAL(8), intent(inout) :: rho(n)
rho(:) = rho(:) + wo*(dble(wf(:))**2 + aimag(wf(:))**2)
RETURN
END SUBROUTINE

END SUBROUTINE 
