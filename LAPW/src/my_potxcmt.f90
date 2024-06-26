SUBROUTINE my_potxcmt(tsh,ias,xctype_,rhomt_,magmt_,taumt_,exmt_,ecmt_,vxcmt_, &
 bxcmt_,wxcmt_)
  ! real input rhomt (density), magmt (magnetization), taumt (kinetic density)
  !
  USE m_atoms, ONLY: natmtot, idxis
  USE m_spin, ONLY: ndmag, nspinor, spinpol, ncmag
  USE m_muffin_tins, ONLY: npmtmax, npmt, nrmt, nrmti
  USE m_density_pot_xc, ONLY: xcgrad, tssxc, ssxc, dncgga, c_tb09
  USE m_states, ONLY: swidth
  USE m_oep_hf, ONLY: hybridc, hybrid
  USE modxcifc, ONLY: xcifc
  USE m_muffin_tins, ONLY: lmmaxo
  !
  IMPLICIT NONE 
  ! arguments
  logical, intent(in) :: tsh
  INTEGER, intent(in) :: ias,xctype_(3)
  REAL(8), intent(in) :: rhomt_(npmtmax,natmtot),magmt_(npmtmax,natmtot,ndmag)
  REAL(8), intent(in) :: taumt_(npmtmax,natmtot,nspinor)
  REAL(8), intent(out) :: exmt_(npmtmax,natmtot),ecmt_(npmtmax,natmtot)
  REAL(8), intent(out) :: vxcmt_(npmtmax,natmtot),bxcmt_(npmtmax,natmtot,ndmag)
  REAL(8), intent(out) :: wxcmt_(npmtmax,natmtot)
  ! local variables
  INTEGER ispn,idm,is
  INTEGER nr,nri,n,i
  REAL(8) t0,t1,t2,t3,t4
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rho(:),rhoup(:),rhodn(:)
  REAL(8), ALLOCATABLE :: gvrho(:,:),gvup(:,:),gvdn(:,:)
  REAL(8), ALLOCATABLE :: grho(:),gup(:),gdn(:)
  REAL(8), ALLOCATABLE :: g2rho(:),g2up(:),g2dn(:)
  REAL(8), ALLOCATABLE :: g3rho(:),g3up(:),g3dn(:)
  REAL(8), ALLOCATABLE :: grho2(:),gup2(:),gdn2(:),gupdn(:)
  REAL(8), ALLOCATABLE :: ex(:),ec(:),vxc(:)
  REAL(8), ALLOCATABLE :: vx(:),vxup(:),vxdn(:)
  REAL(8), ALLOCATABLE :: vc(:),vcup(:),vcdn(:)
  REAL(8), ALLOCATABLE :: mag(:,:),bxc(:,:),tau(:,:)
  REAL(8), ALLOCATABLE :: dxdgr2(:),dxdgu2(:),dxdgd2(:),dxdgud(:)
  REAL(8), ALLOCATABLE :: dcdgr2(:),dcdgu2(:),dcdgd2(:),dcdgud(:)
  REAL(8), ALLOCATABLE :: dxdg2r(:),dxdg2u(:),dxdg2d(:)
  REAL(8), ALLOCATABLE :: dcdg2r(:),dcdg2u(:),dcdg2d(:)
  REAL(8), ALLOCATABLE :: wx(:),wxup(:),wxdn(:)
  REAL(8), ALLOCATABLE :: wc(:),wcup(:),wcdn(:)

  is=idxis(ias)
  n=npmt(is)

  ! allocate local arrays
  ALLOCATE(rho(n),ex(n),ec(n),vxc(n))
  !rho(:) = 0.d0
  !write(*,*) 'rho(1) = ', rho(1)
  !write(*,*) 'rho(2) = ', rho(2)
  !write(*,*) 'sum(rho) = ', sum(rho)
  !stop 'ffr 49'

  IF((xcgrad == 3).or.(xcgrad == 4)) ALLOCATE(tau(n,nspinor))
  IF(spinpol) THEN 
    ALLOCATE(mag(n,3),bxc(n,3))
  ENDIF 
  IF(spinpol) THEN 
    ALLOCATE(rhoup(n),rhodn(n))
    ALLOCATE(vxup(n),vxdn(n),vcup(n),vcdn(n))
    IF(xcgrad == 1) THEN 
      ALLOCATE(grho(n),gup(n),gdn(n))
      ALLOCATE(g2up(n),g2dn(n))
      ALLOCATE(g3rho(n),g3up(n),g3dn(n))
    ELSEIF(xcgrad == 2) THEN 
      ALLOCATE(g2up(n),g2dn(n))
      ALLOCATE(gvup(n,3),gvdn(n,3))
      ALLOCATE(gup2(n),gdn2(n),gupdn(n))
      ALLOCATE(dxdgu2(n),dxdgd2(n),dxdgud(n))
      ALLOCATE(dcdgu2(n),dcdgd2(n),dcdgud(n))
    ELSEIF(xcgrad == 3) THEN 
      ALLOCATE(g2up(n),g2dn(n))
      ALLOCATE(gvup(n,3),gvdn(n,3))
      ALLOCATE(gup2(n),gdn2(n),gupdn(n))
    ELSEIF(xcgrad == 4) THEN 
      ALLOCATE(g2up(n),g2dn(n))
      ALLOCATE(gvup(n,3),gvdn(n,3))
      ALLOCATE(gup2(n),gdn2(n),gupdn(n))
      ALLOCATE(dxdgu2(n),dxdgd2(n),dxdgud(n))
      ALLOCATE(dcdgu2(n),dcdgd2(n),dcdgud(n))
      ALLOCATE(dxdg2u(n),dxdg2d(n))
      ALLOCATE(dcdg2u(n),dcdg2d(n))
      ALLOCATE(wxup(n),wxdn(n),wcup(n),wcdn(n))
    ENDIF 
  ELSE 
    ALLOCATE(vx(n),vc(n))
    IF(xcgrad == 1) THEN 
      ALLOCATE(grho(n),g2rho(n),g3rho(n))
    ELSEIF(xcgrad == 2) THEN 
      ALLOCATE(g2rho(n),gvrho(n,3),grho2(n))
      ALLOCATE(dxdgr2(n),dcdgr2(n))
    ELSEIF(xcgrad == 3) THEN 
      ALLOCATE(g2rho(n),gvrho(n,3),grho2(n))
    ELSEIF(xcgrad == 4) THEN 
      ALLOCATE(g2rho(n),gvrho(n,3),grho2(n))
      ALLOCATE(dxdgr2(n),dcdgr2(n))
      ALLOCATE(dxdg2r(n),dcdg2r(n))
      ALLOCATE(wx(n),wc(n))
    ENDIF 
  ENDIF 

  write(*,*)
  write(*,*) 'my_potxcmt: xcgrad = ', xcgrad
  write(*,*) 'my_potxcmt: xctype_ = ', xctype_
  write(*,*) 'my_potxcmt: shape(rho) = ', shape(rho)
  write(*,*) 'my_potxcmt: shape(rhomt) = ', shape(rhomt_)
  write(*,'(1x,A,ES18.10)') 'sum(rho)    = ', sum(rho)
  write(*,'(1x,A,ES18.10)') 'sum(rhomt_) = ', sum(rhomt_)

  !
  nr=nrmt(is)
  nri=nrmti(is)
  !
  IF(tsh) THEN 
    ! convert the density to spherical coordinates
    CALL my_rbsht(nr,nri,rhomt_(:,ias),rho)
  ELSE 
    rho(1:n)=rhomt_(1:n,ias)
  ENDIF 

  !write(*,'(1x,A,ES18.10)') 'sum(rho)    = ', sum(rho)
  !write(*,'(1x,A,ES18.10)') 'sum(rhomt_) = ', sum(rhomt_)

  !is = 1
  !write(*,*) 'Some rhomt after my_rbsht'
  !write(*,*) 'Inner'
  !do i = 1,10
  !  write(*,'(1x,I4,ES18.10)') i, rho(i)
  !enddo
  !write(*,*) 'Outer'
  !do i = nrmti(is)+1,nrmti(is)+lmmaxo
  !  write(*,'(1x,I4,ES18.10)') i, rho(i)
  !enddo


  ! convert tau to spherical coordinates if required
  IF((xcgrad == 3).or.(xcgrad == 4)) THEN 
    DO ispn=1,nspinor
      IF(tsh) THEN 
        CALL my_rbsht(nr,nri,taumt_(:,ias,ispn),tau(:,ispn))
      ELSE 
        tau(1:n,ispn)=taumt_(1:n,ias,ispn)
      ENDIF 
    ENDDO 
  ENDIF 

  IF(spinpol) THEN 
  !------------------------!
  !     spin-polarised     !
  !------------------------!
  ! magnetisation in spherical coordinates
    DO idm=1,ndmag
      IF(tsh) THEN 
        CALL my_rbsht(nr,nri,magmt_(:,ias,idm),mag(:,idm))
      ELSE 
        mag(1:n,idm)=magmt_(1:n,ias,idm)
      ENDIF 
    ENDDO 
    ! use scaled spin exchange-correlation (SSXC) if required
    IF(tssxc) mag(:,1:ndmag)=mag(:,1:ndmag)*ssxc
    !
    IF(ncmag) THEN 
      ! non-collinear (use Kubler's trick)
      IF(xcgrad == 0) THEN 
        ! LSDA
        DO i=1,n
          ! compute rhoup=(rho+|m|)/2 and rhodn=(rho-|m|)/2
          t0=rho(i)
          t1=sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2)
          rhoup(i)=0.5d0*(t0+t1)
          rhodn(i)=0.5d0*(t0-t1)
        ENDDO 
      ELSE 
        ! functionals which require gradients
        DO i=1,n
          t0=rho(i)
          t1=sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2+dncgga)
          rhoup(i)=0.5d0*(t0+t1)
          rhodn(i)=0.5d0*(t0-t1)
        ENDDO 
      ENDIF 
    ELSE 
      ! collinear
      DO i=1,n
        ! compute rhoup=(rho+m_z)/2 and rhodn=(rho-m_z)/2
        t0=rho(i)
        t1=mag(i,1)
        rhoup(i)=0.5d0*(t0+t1)
        rhodn(i)=0.5d0*(t0-t1)
      ENDDO 
    ENDIF 
    ! CALL the exchange-correlation interface routine
    IF(xcgrad <= 0) THEN 
      CALL xcifc(xctype_,n=n,tempa=swidth,rhoup=rhoup,rhodn=rhodn,ex=ex,ec=ec, &
       vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
    !
    ELSEIF(xcgrad == 1) THEN 
      CALL ggamt_sp_1(is,n,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
      CALL xcifc(xctype_,n=n,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup,gdn=gdn, &
       g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=ex,ec=ec,vxup=vxup,&
       vxdn=vxdn,vcup=vcup,vcdn=vcdn)
    !
    ELSEIF(xcgrad == 2) THEN 
      CALL ggamt_sp_2a(is,n,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
      CALL xcifc(xctype_,n=n,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2, &
       gupdn=gupdn,ex=ex,ec=ec,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn, &
       dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2, &
       dcdgud=dcdgud)
      CALL ggamt_sp_2b(is,n,g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2, &
       dxdgd2,dxdgud,dcdgu2,dcdgd2,dcdgud)
    !
    ELSEIF(xcgrad == 3) THEN 
      CALL ggamt_sp_2a(is,n,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
      CALL xcifc(xctype_,n=n,c_tb09=c_tb09,rhoup=rhoup,rhodn=rhodn,g2up=g2up, &
       g2dn=g2dn,gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=tau(:,1),taudn=tau(:,2), &
       vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
       ex(:)=0.d0; ec(:)=0.d0
    !
    ELSEIF(xcgrad == 4) THEN 
      CALL ggamt_sp_2a(is,n,rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
      CALL xcifc(xctype_,n=n,rhoup=rhoup,rhodn=rhodn,g2up=g2up,g2dn=g2dn, &
       gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=tau(:,1),taudn=tau(:,2),ex=ex,ec=ec,&
       vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn,dxdgu2=dxdgu2,dxdgd2=dxdgd2, &
       dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2,dcdgud=dcdgud,dxdg2u=dxdg2u, &
       dxdg2d=dxdg2d,dcdg2u=dcdg2u,dcdg2d=dcdg2d,wxup=wxup,wxdn=wxdn,wcup=wcup, &
       wcdn=wcdn)
      CALL ggamt_sp_2b(is,n,g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2, &
       dxdgd2,dxdgud,dcdgu2,dcdgd2,dcdgud)
      wxup(:)=0.5d0*(wxup(:)+wxdn(:)+wcup(:)+wcdn(:))
      IF(tsh) THEN 
        CALL my_rfsht(nr,nri,wxup,wxcmt_(:,ias))
      ELSE 
        wxcmt_(1:n,ias)=wxup(1:n)
      ENDIF 
    ENDIF 
    ! hybrid functionals
    IF(hybrid) THEN 
      t1=1.d0-hybridc
      ! scale exchange part of energy
      ex(:)=t1*ex(:)
      ! scale exchange part of potential
      vxup(:)=t1*vxup(:)
      vxdn(:)=t1*vxdn(:)
    ENDIF 
    IF(ncmag) THEN 
      ! non-collinear: locally spin rotate the exchange-correlation potential
      DO i=1,n
        t1=vxup(i)+vcup(i)
        t2=vxdn(i)+vcdn(i)
        vxc(i)=0.5d0*(t1+t2)
        ! determine the exchange-correlation magnetic field
        t3=0.5d0*(t1-t2)
        ! |m| = rhoup - rhodn
        t4=rhoup(i)-rhodn(i)
        IF(abs(t4).gt.1.d-8) t4=t3/t4
        bxc(i,1:3)=mag(i,1:3)*t4
      ENDDO 
    ELSE 
      ! collinear
      DO i=1,n
        t1=vxup(i)+vcup(i)
        t2=vxdn(i)+vcdn(i)
        vxc(i)=0.5d0*(t1+t2)
        bxc(i,1)=0.5d0*(t1-t2)
      ENDDO 
    ENDIF 
    ! scale B_xc for SSXC if required
    IF(tssxc) bxc(:,1:ndmag)=bxc(:,1:ndmag)*ssxc
    DO idm=1,ndmag
      IF(tsh) THEN 
        ! convert field to spherical harmonics
        CALL my_rfsht(nr,nri,bxc(:,idm),bxcmt_(:,ias,idm))
      ELSE 
        bxcmt_(1:n,ias,idm)=bxc(1:n,idm)
      ENDIF 
    ENDDO 
  
  ELSE 

    !--------------------------!
    !     spin-unpolarised     !
    !--------------------------!

    IF(xcgrad <= 0) THEN 
      CALL xcifc(xctype_,n=n,tempa=swidth,rho=rho,ex=ex,ec=ec,vx=vx,vc=vc)
    ELSEIF(xcgrad == 1) THEN 
      CALL ggamt_1(tsh,is,n,rhomt_(:,ias),grho,g2rho,g3rho)
      CALL xcifc(xctype_,n=n,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, &
       ec=ec,vx=vx,vc=vc)
    ELSEIF(xcgrad == 2) THEN 
      CALL ggamt_2a(tsh,is,n,rhomt_(:,ias),g2rho,gvrho,grho2)
      CALL xcifc(xctype_,n=n,rho=rho,grho2=grho2,ex=ex,ec=ec,vx=vx,vc=vc, &
       dxdgr2=dxdgr2,dcdgr2=dcdgr2)
      CALL ggamt_2b(is,n,g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
    ELSEIF(xcgrad == 3) THEN 
      CALL ggamt_2a(tsh,is,n,rhomt_(:,ias),g2rho,gvrho,grho2)
      CALL xcifc(xctype_,n=n,c_tb09=c_tb09,rho=rho,g2rho=g2rho,grho2=grho2, &
       tau=tau,vx=vx,vc=vc)
      ex(:)=0.d0; ec(:)=0.d0
    ELSEIF(xcgrad == 4) THEN 
      CALL ggamt_2a(tsh,is,n,rhomt_(:,ias),g2rho,gvrho,grho2)
      CALL xcifc(xctype_,n=n,rho=rho,g2rho=g2rho,grho2=grho2,tau=tau,ex=ex,ec=ec,&
       vx=vx,vc=vc,dxdgr2=dxdgr2,dcdgr2=dcdgr2,dxdg2r=dxdg2r,dcdg2r=dcdg2r,wx=wx,&
       wc=wc)
      CALL ggamt_2b(is,n,g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
      wx(:)=wx(:)+wc(:)
      IF(tsh) THEN 
        CALL my_rfsht(nr,nri,wx,wxcmt_(:,ias))
      ELSE 
        wxcmt_(1:n,ias)=wx(1:n)
      ENDIF 
    ENDIF 
    ! hybrid functionals
    IF(hybrid) THEN 
      t1=1.d0-hybridc
      ! scale exchange part of energy
      ex(:)=t1*ex(:)
      ! scale exchange part of potential
      vxc(:)=t1*vx(:)+vc(:)
    ELSE 
      vxc(:)=vx(:)+vc(:)
    ENDIF 
  ENDIF
  write(*,*)
  write(*,*) 'my_potxcmt: tsh = ', tsh 
  IF(tsh) THEN 
    ! convert exchange and correlation energy densities to spherical harmonics
    CALL my_rfsht(nr,nri,ex,exmt_(:,ias))
    CALL my_rfsht(nr,nri,ec,ecmt_(:,ias))
    ! convert exchange-correlation potential to spherical harmonics
    CALL my_rfsht(nr,nri,vxc,vxcmt_(:,ias))

    write(*,*) 'my_potxcmt: shape(exmt)  = ', shape(exmt_)
    write(*,*) 'my_potxcmt: shape(ecmt)  = ', shape(ecmt_)
    write(*,*) 'my_potxcmt: shape(vxcmt) = ', shape(vxcmt_)

  ELSE 
    exmt_(1:n,ias)=ex(1:n)
    ecmt_(1:n,ias)=ec(1:n)
    vxcmt_(1:n,ias)=vxc(1:n)
  ENDIF 

  ! Deallocate memory
  DEALLOCATE(rho,ex,ec,vxc)
  IF((xcgrad == 3).or.(xcgrad == 4)) DEALLOCATE(tau)
  IF(spinpol) THEN 
    DEALLOCATE(mag,bxc)
    DEALLOCATE(rhoup,rhodn,vxup,vxdn,vcup,vcdn)
    IF(xcgrad == 1) THEN 
      DEALLOCATE(grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
    ELSEIF(xcgrad == 2) THEN 
      DEALLOCATE(g2up,g2dn)
      DEALLOCATE(gvup,gvdn)
      DEALLOCATE(gup2,gdn2,gupdn)
      DEALLOCATE(dxdgu2,dxdgd2,dxdgud)
      DEALLOCATE(dcdgu2,dcdgd2,dcdgud)
    ELSEIF(xcgrad == 3) THEN 
      DEALLOCATE(g2up,g2dn)
      DEALLOCATE(gvup,gvdn)
      DEALLOCATE(gup2,gdn2,gupdn)
    ENDIF 
  ELSE 
    DEALLOCATE(vx,vc)
    IF(xcgrad == 1) THEN 
      DEALLOCATE(grho,g2rho,g3rho)
    ELSEIF(xcgrad == 2) THEN 
      DEALLOCATE(g2rho,gvrho,grho2)
      DEALLOCATE(dxdgr2,dcdgr2)
    ELSEIF(xcgrad == 3) THEN 
      DEALLOCATE(g2rho,gvrho,grho2)
    ELSEIF(xcgrad == 4) THEN 
      DEALLOCATE(g2rho,gvrho,grho2)
      DEALLOCATE(dxdgr2,dcdgr2,dxdg2r,dcdg2r)
      DEALLOCATE(wx,wc)
    ENDIF 
  ENDIF 
  RETURN 
END SUBROUTINE 

