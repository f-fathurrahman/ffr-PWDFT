SUBROUTINE my_potxcir(xctype_,rhoir_,magir_,tauir_,exir_,ecir_,vxcir_,bxcir_,wxcir_)
  !
  USE m_gvectors, ONLY: ngtot
  USE m_spin, ONLY: ndmag, nspinor, spinpol, ncmag
  USE m_density_pot_xc, ONLY: xcgrad, tssxc, ssxc, dncgga, c_tb09
  USE m_states, ONLY: swidth
  USE m_oep_hf, ONLY: hybridc, hybrid
  USE modxcifc, ONLY: xcifc
  !
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: xctype_(3)
  REAL(8), intent(in) :: rhoir_(ngtot),magir_(ngtot,ndmag),tauir_(ngtot,nspinor)
  REAL(8), intent(out) :: exir_(ngtot),ecir_(ngtot)
  REAL(8), intent(out) :: vxcir_(ngtot),bxcir_(ngtot,ndmag),wxcir_(ngtot)
  ! local variables
  INTEGER n,i
  REAL(8) t0,t1,t2,t3,t4
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rhoup(:),rhodn(:)
  REAL(8), ALLOCATABLE :: gvrho(:,:),gvup(:,:),gvdn(:,:)
  REAL(8), ALLOCATABLE :: grho(:),gup(:),gdn(:)
  REAL(8), ALLOCATABLE :: g2rho(:),g2up(:),g2dn(:)
  REAL(8), ALLOCATABLE :: g3rho(:),g3up(:),g3dn(:)
  REAL(8), ALLOCATABLE :: grho2(:),gup2(:),gdn2(:),gupdn(:)
  REAL(8), ALLOCATABLE :: vx(:),vxup(:),vxdn(:)
  REAL(8), ALLOCATABLE :: vc(:),vcup(:),vcdn(:)
  REAL(8), ALLOCATABLE :: dxdgr2(:),dxdgu2(:),dxdgd2(:),dxdgud(:)
  REAL(8), ALLOCATABLE :: dcdgr2(:),dcdgu2(:),dcdgd2(:),dcdgud(:)
  REAL(8), ALLOCATABLE :: dxdg2r(:),dxdg2u(:),dxdg2d(:)
  REAL(8), ALLOCATABLE :: dcdg2r(:),dcdg2u(:),dcdg2d(:)
  REAL(8), ALLOCATABLE :: wx(:),wxup(:),wxdn(:)
  REAL(8), ALLOCATABLE :: wc(:),wcup(:),wcdn(:)
  n=ngtot
  IF(spinpol) THEN 
    ALLOCATE(rhoup(n),rhodn(n))
    ALLOCATE(vxup(n),vxdn(n),vcup(n),vcdn(n))
    IF(xcgrad.eq.1) THEN 
      ALLOCATE(grho(n),gup(n),gdn(n))
      ALLOCATE(g2up(n),g2dn(n))
      ALLOCATE(g3rho(n),g3up(n),g3dn(n))
    ELSEIF(xcgrad.eq.2) THEN 
      ALLOCATE(g2up(n),g2dn(n))
      ALLOCATE(gvup(n,3),gvdn(n,3))
      ALLOCATE(gup2(n),gdn2(n),gupdn(n))
      ALLOCATE(dxdgu2(n),dxdgd2(n),dxdgud(n))
      ALLOCATE(dcdgu2(n),dcdgd2(n),dcdgud(n))
    ELSEIF(xcgrad.eq.3) THEN 
      ALLOCATE(g2up(n),g2dn(n))
      ALLOCATE(gvup(n,3),gvdn(n,3))
      ALLOCATE(gup2(n),gdn2(n),gupdn(n))
    ELSEIF(xcgrad.eq.4) THEN 
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
    IF(xcgrad.eq.1) THEN 
      ALLOCATE(grho(n),g2rho(n),g3rho(n))
    ELSEIF(xcgrad.eq.2) THEN 
      ALLOCATE(g2rho(n),gvrho(n,3),grho2(n))
      ALLOCATE(dxdgr2(n),dcdgr2(n))
    ELSEIF(xcgrad.eq.3) THEN 
      ALLOCATE(g2rho(n),gvrho(n,3),grho2(n))
    ELSEIF(xcgrad.eq.4) THEN 
      ALLOCATE(g2rho(n),gvrho(n,3),grho2(n))
      ALLOCATE(dxdgr2(n),dcdgr2(n))
      ALLOCATE(dxdg2r(n),dcdg2r(n))
      ALLOCATE(wx(n),wc(n))
    ENDIF 
  ENDIF 
  IF(spinpol) THEN 
  !------------------------!
  !     spin-polarised     !
  !------------------------!
    IF(ncmag) THEN 
  ! non-collinear
      IF(xcgrad.eq.0) THEN 
  ! LSDA
        DO i=1,n
          t0=rhoir_(i)
          t1=sqrt(magir_(i,1)**2+magir_(i,2)**2+magir_(i,3)**2)*ssxc
          rhoup(i)=0.5d0*(t0+t1)
          rhodn(i)=0.5d0*(t0-t1)
        ENDDO 
      ELSE 
  ! functionals which require gradients
        DO i=1,n
          t0=rhoir_(i)
          t1=sqrt(magir_(i,1)**2+magir_(i,2)**2+magir_(i,3)**2+dncgga)*ssxc
          rhoup(i)=0.5d0*(t0+t1)
          rhodn(i)=0.5d0*(t0-t1)
        ENDDO 
      ENDIF 
    ELSE 
  ! collinear
      DO i=1,n
        t0=rhoir_(i)
        t1=magir_(i,1)*ssxc
        rhoup(i)=0.5d0*(t0+t1)
        rhodn(i)=0.5d0*(t0-t1)
      ENDDO 
    ENDIF 
    IF(xcgrad.le.0) THEN 
      CALL xcifc(xctype_,n=n,tempa=swidth,rhoup=rhoup,rhodn=rhodn,ex=exir_, &
       ec=ecir_,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
    ELSEIF(xcgrad.eq.1) THEN 
      CALL ggair_sp_1(rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
      CALL xcifc(xctype_,n=n,rhoup=rhoup,rhodn=rhodn,grho=grho,gup=gup,gdn=gdn, &
       g2up=g2up,g2dn=g2dn,g3rho=g3rho,g3up=g3up,g3dn=g3dn,ex=exir_,ec=ecir_, &
       vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
    ELSEIF(xcgrad.eq.2) THEN 
      CALL ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
      CALL xcifc(xctype_,n=n,rhoup=rhoup,rhodn=rhodn,gup2=gup2,gdn2=gdn2, &
       gupdn=gupdn,ex=exir_,ec=ecir_,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn, &
       dxdgu2=dxdgu2,dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2, &
       dcdgud=dcdgud)
      CALL ggair_sp_2b(g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2,dxdgd2, &
       dxdgud,dcdgu2,dcdgd2,dcdgud)
    ELSEIF(xcgrad.eq.3) THEN 
      CALL ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
      CALL xcifc(xctype_,n=n,c_tb09=c_tb09,rhoup=rhoup,rhodn=rhodn,g2up=g2up, &
       g2dn=g2dn,gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=tauir_(:,1), &
       taudn=tauir_(:,2),vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn)
      exir_(:)=0.d0; ecir_(:)=0.d0
    ELSEIF(xcgrad.eq.4) THEN 
      CALL ggair_sp_2a(rhoup,rhodn,g2up,g2dn,gvup,gvdn,gup2,gdn2,gupdn)
      CALL xcifc(xctype_,n=n,rhoup=rhoup,rhodn=rhodn,g2up=g2up,g2dn=g2dn, &
       gup2=gup2,gdn2=gdn2,gupdn=gupdn,tauup=tauir_(:,1),taudn=tauir_(:,2), &
       ex=exir_,ec=ecir_,vxup=vxup,vxdn=vxdn,vcup=vcup,vcdn=vcdn,dxdgu2=dxdgu2, &
       dxdgd2=dxdgd2,dxdgud=dxdgud,dcdgu2=dcdgu2,dcdgd2=dcdgd2,dcdgud=dcdgud, &
       dxdg2u=dxdg2u,dxdg2d=dxdg2d,dcdg2u=dcdg2u,dcdg2d=dcdg2d,wxup=wxup, &
       wxdn=wxdn,wcup=wcup,wcdn=wcdn)
      CALL ggair_sp_2b(g2up,g2dn,gvup,gvdn,vxup,vxdn,vcup,vcdn,dxdgu2,dxdgd2, &
       dxdgud,dcdgu2,dcdgd2,dcdgud)
      wxcir_(:)=0.5d0*(wxup(:)+wxdn(:)+wcup(:)+wcdn(:))
    ENDIF 
  ! hybrid functionals
    IF(hybrid) THEN 
      t1=1.d0-hybridc
  ! scale exchange part of energy
      exir_(:)=t1*exir_(:)
  ! scale exchange part of potential
      vxup(:)=t1*vxup(:)
      vxdn(:)=t1*vxdn(:)
    ENDIF 
    IF(ncmag) THEN 
  ! non-collinear: spin rotate the local exchange potential
      DO i=1,n
        t1=vxup(i)+vcup(i)
        t2=vxdn(i)+vcdn(i)
        vxcir_(i)=0.5d0*(t1+t2)
  ! determine the exchange-correlation magnetic field
        t3=0.5d0*(t1-t2)
        t4=rhoup(i)-rhodn(i)
        IF(abs(t4).gt.1.d-8) t4=t3/t4
        bxcir_(i,:)=magir_(i,:)*t4
      ENDDO 
    ELSE 
  ! collinear
      DO i=1,n
        t1=vxup(i)+vcup(i)
        t2=vxdn(i)+vcdn(i)
        vxcir_(i)=0.5d0*(t1+t2)
        bxcir_(i,1)=0.5d0*(t1-t2)
      ENDDO 
    ENDIF 
  ! scale field if required
    IF(tssxc) bxcir_(:,1:ndmag)=bxcir_(:,1:ndmag)*ssxc
  ELSE 
  !--------------------------!
  !     spin-unpolarised     !
  !--------------------------!
    IF(xcgrad.le.0) THEN 
      CALL xcifc(xctype_,n=n,tempa=swidth,rho=rhoir_,ex=exir_,ec=ecir_,vx=vx, &
       vc=vc)
    ELSEIF(xcgrad.eq.1) THEN 
      CALL ggair_1(rhoir_,grho,g2rho,g3rho)
      CALL xcifc(xctype_,n=n,rho=rhoir_,grho=grho,g2rho=g2rho,g3rho=g3rho, &
       ex=exir_,ec=ecir_,vx=vx,vc=vc)
    ELSEIF(xcgrad.eq.2) THEN 
      CALL ggair_2a(rhoir_,g2rho,gvrho,grho2)
      CALL xcifc(xctype_,n=n,rho=rhoir_,grho2=grho2,ex=exir_,ec=ecir_,vx=vx, &
       vc=vc,dxdgr2=dxdgr2,dcdgr2=dcdgr2)
      CALL ggair_2b(g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
    ELSEIF(xcgrad.eq.3) THEN 
      CALL ggair_2a(rhoir_,g2rho,gvrho,grho2)
      CALL xcifc(xctype_,n=n,c_tb09=c_tb09,rho=rhoir_,g2rho=g2rho,grho2=grho2, &
       tau=tauir_,vx=vx,vc=vc)
      exir_(:)=0.d0; ecir_(:)=0.d0
    ELSEIF(xcgrad.eq.4) THEN 
      CALL ggair_2a(rhoir_,g2rho,gvrho,grho2)
      CALL xcifc(xctype_,n=n,rho=rhoir_,g2rho=g2rho,grho2=grho2,tau=tauir_, &
       ex=exir_,ec=ecir_,vx=vx,vc=vc,dxdgr2=dxdgr2,dcdgr2=dcdgr2,dxdg2r=dxdg2r, &
       dcdg2r=dcdg2r,wx=wx,wc=wc)
      CALL ggair_2b(g2rho,gvrho,vx,vc,dxdgr2,dcdgr2)
      wxcir_(:)=wx(:)+wc(:)
    ENDIF 
  ! hybrid functionals
    IF(hybrid) THEN 
      t1=1.d0-hybridc
  ! scale exchange part of energy
      exir_(:)=t1*exir_(:)
  ! scale exchange part of potential
      vxcir_(:)=t1*vx(:)+vc(:)
    ELSE 
      vxcir_(:)=vx(:)+vc(:)
    ENDIF 
  ENDIF 
  IF(spinpol) THEN 
    DEALLOCATE(rhoup,rhodn,vxup,vxdn,vcup,vcdn)
    IF(xcgrad.eq.1) THEN 
      DEALLOCATE(grho,gup,gdn,g2up,g2dn,g3rho,g3up,g3dn)
    ELSEIF(xcgrad.eq.2) THEN 
      DEALLOCATE(g2up,g2dn)
      DEALLOCATE(gvup,gvdn)
      DEALLOCATE(gup2,gdn2,gupdn)
      DEALLOCATE(dxdgu2,dxdgd2,dxdgud)
      DEALLOCATE(dcdgu2,dcdgd2,dcdgud)
    ELSEIF(xcgrad.eq.3) THEN 
      DEALLOCATE(g2up,g2dn)
      DEALLOCATE(gvup,gvdn)
      DEALLOCATE(gup2,gdn2,gupdn)
    ENDIF 
  ELSE 
    DEALLOCATE(vx,vc)
    IF(xcgrad.eq.1) THEN 
      DEALLOCATE(grho,g2rho,g3rho)
    ELSEIF(xcgrad.eq.2) THEN 
      DEALLOCATE(g2rho,gvrho,grho2)
      DEALLOCATE(dxdgr2,dcdgr2)
    ELSEIF(xcgrad.eq.3) THEN 
      DEALLOCATE(g2rho,gvrho,grho2)
    ELSEIF(xcgrad.eq.4) THEN 
      DEALLOCATE(g2rho,gvrho,grho2)
      DEALLOCATE(dxdgr2,dcdgr2,dxdg2r,dcdg2r)
      DEALLOCATE(wx,wc)
    ENDIF 
  ENDIF 
  RETURN 
END SUBROUTINE 

