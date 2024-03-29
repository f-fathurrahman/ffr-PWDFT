SUBROUTINE rhoinit()
  USE m_density_pot_xc, ONLY: rhoir, rhomt
  USE m_muffin_tins, ONLY: nrmti, nrmt, rcmt, rcmt, npcmtmax, npcmt, nrcmt, nrcmti, nrcmtmax, &
                   lmmaxi, lmmaxo, lmaxi
  USE m_constants, ONLY: zil, y00, fourpi
  USE m_lattice, ONLY: omega, epslat
  USE m_atoms, ONLY: idxis, idxas, nspecies, natoms, natmtot
  USE m_atomic_species, ONLY: nrsp, nrspmax, rsp, nrsp, rhosp
  USE m_gvectors, ONLY: igfft, sfacg, gc, gmaxvr, ngvec, ngtot, ngridg, ylmg
  USE m_charge_moment_current, ONLY: chgexs
  IMPLICIT NONE 
  ! local variables
  INTEGER :: lmax,is,ia,ias,i
  INTEGER :: nr,nri,nro,nrs,ir
  INTEGER :: nrc,nrci,irco,irc
  INTEGER :: l,m,lm,ig,ifg
  REAL(8) :: x,t1,t2
  COMPLEX(8) :: z1,z2,z3
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: jl(:,:),ffg(:),wr(:),fr(:)
  COMPLEX(8), ALLOCATABLE :: zfmt(:),zfft(:)
  
  lmax=min(lmaxi,1)
  ! zero the charge density arrays
  rhomt(:,:)=0.d0
  rhoir(:)=0.d0
  ! compute the superposition of all the atomic density tails
  ALLOCATE(zfft(ngtot))
  zfft(:)=0.d0
  
  ALLOCATE(ffg(ngvec),wr(nrspmax),fr(nrspmax))

  DO is=1,nspecies
    nr=nrmt(is)
    nrs=nrsp(is)
    nro=nrs-nr+1
    ! determine the weights for the radial integral
    CALL wsplint(nro,rsp(nr,is),wr(nr))
    DO ig=1,ngvec
      t1=gc(ig)
      ! spherical bessel function j_0(x) times the atomic density tail
      IF(t1.gt.epslat) THEN 
        t2=1.d0/t1
        DO ir=nr,nrs
          x=t1*rsp(ir,is)
          fr(ir)=t2*sin(x)*rhosp(ir,is)*rsp(ir,is)
        ENDDO 
      ELSE 
        fr(nr:nrs)=rhosp(nr:nrs,is)*rsp(nr:nrs,is)**2
      ENDIF 
      t1=dot_product(wr(nr:nrs),fr(nr:nrs))
      ! apply low-pass filter
      t1=t1*exp(-4.d0*(gc(ig)/gmaxvr)**2)
      ffg(ig)=(fourpi/omega)*t1
    ENDDO 
    DO ia=1,natoms(is)
      ias=idxas(ia,is)
      DO ig=1,ngvec
        ifg=igfft(ig)
        zfft(ifg)=zfft(ifg)+ffg(ig)*conjg(sfacg(ig,ias))
      ENDDO 
    ENDDO 
  ENDDO 
  DEALLOCATE(ffg,wr,fr)

  ! compute the tails in each muffin-tin
  ALLOCATE(jl(0:lmax,nrcmtmax),zfmt(npcmtmax))
  DO ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    irco=nrci+1
    zfmt(1:npcmt(is))=0.d0
    DO ig=1,ngvec
      ifg=igfft(ig)
      DO irc=1,nrc
        x=gc(ig)*rcmt(irc,is)
        CALL sbessel(lmax,x,jl(:,irc))
      ENDDO 
      z1=fourpi*zfft(ifg)*sfacg(ig,ias)
      lm=0
      DO l=0,lmax
        z2=z1*zil(l)
        DO m=-l,l
          lm=lm+1
          z3=z2*conjg(ylmg(lm,ig))
          i=lm
          DO irc=1,nrci
            zfmt(i)=zfmt(i)+jl(l,irc)*z3
            i=i+lmmaxi
          ENDDO 
          DO irc=irco,nrc
            zfmt(i)=zfmt(i)+jl(l,irc)*z3
            i=i+lmmaxo
          ENDDO 
        ENDDO 
      ENDDO 
    ENDDO 
    CALL ztorfmt(nrc,nrci,zfmt,rhomt(:,ias))
  ENDDO 
  DEALLOCATE(jl,zfmt)

  ! convert the density from a coarse to a fine radial mesh
  CALL rfmtctof(rhomt)

  ! add the atomic charge density and the excess charge in each muffin-tin
  t1=chgexs/omega
  DO ias=1,natmtot
    is=idxis(ias)
    nr=nrmt(is)
    nri=nrmti(is)
    i=1
    DO ir=1,nri
      t2=(t1+rhosp(ir,is))/y00
      rhomt(i,ias)=rhomt(i,ias)+t2
      i=i+lmmaxi
    ENDDO 
    DO ir=nri+1,nr
      t2=(t1+rhosp(ir,is))/y00
      rhomt(i,ias)=rhomt(i,ias)+t2
      i=i+lmmaxo
    ENDDO 
  ENDDO 

  ! interstitial density determined from the atomic tails and excess charge
  CALL zfftifc(3,ngridg,1,zfft)
  DO ir=1,ngtot
    rhoir(ir)=dble(zfft(ir))+t1
    ! make sure that the density is always positive
    IF(rhoir(ir).lt.1.d-10) rhoir(ir)=1.d-10
  ENDDO 
  DEALLOCATE(zfft)
    
  RETURN 
END SUBROUTINE 