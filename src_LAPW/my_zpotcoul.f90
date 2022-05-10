SUBROUTINE my_zpotcoul(nr,nri,np,npi,ld1,rl,ngdg,igf,ngp,gpc,gclgp,ld2,jlgprmt, &
  ylmgp,sfacgp,zrhoir,ld3,zvclmt,zvclir)
  !
  USE m_constants, ONLY: zil, zilc, y00, fourpi
  USE m_atoms, ONLY: idxas, idxis, natoms, natmtot, nspecies
  USE m_muffin_tins, ONLY: lmaxi, lmmaxo, rmt, lmaxo, lmmaxi
  USE m_lattice, ONLY: omega, epslat
  USE m_density_pot_xc, ONLY: lnpsd
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr(nspecies),nri(nspecies)
  INTEGER, intent(in) :: np(nspecies),npi(nspecies)
  INTEGER, intent(in) :: ld1
  REAL(8), intent(in) :: rl(ld1,-lmaxo-1:lmaxo+2,nspecies)
  INTEGER, intent(in) :: ngdg(3),igf(*),ngp
  REAL(8), intent(in) :: gpc(ngp),gclgp(ngp)
  INTEGER, intent(in) :: ld2
  REAL(8), intent(in) :: jlgprmt(0:lnpsd,ld2,nspecies)
  COMPLEX(8), intent(in) :: ylmgp(lmmaxo,ngp),sfacgp(ld2,natmtot)
  COMPLEX(8), intent(in) :: zrhoir(*)
  INTEGER, intent(in) :: ld3
  COMPLEX(8), intent(inout) :: zvclmt(ld3,natmtot)
  COMPLEX(8), intent(out) :: zvclir(*)
  ! local variables
  INTEGER is,ia,ias
  INTEGER ir,l,m,lm
  INTEGER ig,jg,i,j
  REAL(8) t1,t2,t3
  COMPLEX(8) z1,z2
  ! automatic arrays
  REAL(8) :: rmtl(0:lmaxo+3,nspecies)
  COMPLEX(8) :: qlm(lmmaxo,natmtot)
  COMPLEX(8) :: zl(0:lmaxo),zlm(lmmaxo)
  COMPLEX(8) :: zhmt(ld3)
  ! external functions
  REAL(8) :: factnm
  external factnm

  ! compute (R_mt)^l
  DO is = 1,nspecies
    rmtl(0,is)=1.d0
    DO l = 1,lmaxo+3
      rmtl(l,is) = rmtl(l-1,is)*rmt(is)
    ENDDO 
  ENDDO 

  ! compute the multipole moments from the muffin-tin potentials
  t1 = 1.d0/fourpi
  DO ias = 1,natmtot
    is = idxis(ias)
    i = np(is) - lmmaxo
    lm = 0
    DO l = 0,lmaxo
      t2 = t1*dble(2*l+1)*rmtl(l+1,is)
      DO m = -l,l
        lm = lm + 1
        i = i + 1
        qlm(lm,ias) = t2*zvclmt(i,ias)
      ENDDO 
    ENDDO 
  ENDDO 

  ! Fourier transform density to G-space and store in zvclir
  CALL zcopy(ngdg(1)*ngdg(2)*ngdg(3),zrhoir,1,zvclir,1)
  CALL zfftifc(3,ngdg,-1,zvclir)

  ! subtract the multipole moments of the interstitial charge density
  DO is = 1,nspecies
    DO l = 0,lmaxo
      zl(l) = fourpi*zil(l)*rmtl(l+2,is)
    ENDDO 
    DO ia = 1,natoms(is)
      ias = idxas(ia,is)
      DO ig = 1,ngp
        jg = igf(ig)
        IF(gpc(ig) .gt. epslat) THEN 
          z1=zvclir(jg)*sfacgp(ig,ias)/gpc(ig)
          lm=0
          DO l=0,lmaxo
            z2=jlgprmt(l+1,ig,is)*z1*zl(l)
            DO m=-l,l
              lm=lm+1
              qlm(lm,ias)=qlm(lm,ias)-z2*conjg(ylmgp(lm,ig))
            ENDDO 
          ENDDO 
        ELSE 
          t1=(fourpi/3.d0)*rmtl(3,is)*y00
          qlm(1,ias)=qlm(1,ias)-t1*zvclir(jg)
        ENDIF 
      ENDDO 
    ENDDO 
  ENDDO 

  ! find the smooth pseudocharge within the muffin-tin whose multipoles are the
  ! difference between the real muffin-tin and interstitial multipoles
  t1=(fourpi/omega)*factnm(2*lnpsd+1,2)
  DO ias=1,natmtot
    is=idxis(ias)
    lm=0
    DO l=0,lmaxo
      t2=t1/(factnm(2*l+1,2)*rmtl(l,is))
      z1=t2*zilc(l)
      DO m=-l,l
        lm=lm+1
        zlm(lm)=z1*qlm(lm,ias)
      ENDDO 
    ENDDO 
    
    ! add the pseudocharge and real interstitial densities in G-space
    DO ig=1,ngp
      jg=igf(ig)
      IF(gpc(ig).gt.epslat) THEN 
        t2=gpc(ig)*rmt(is)
        t3=1.d0/t2**lnpsd
        z1=t3*zlm(1)*ylmgp(1,ig)
        lm=1
        DO l=1,lmaxo
          lm=lm+1
          z2=zlm(lm)*ylmgp(lm,ig)
          DO m=1-l,l
            lm=lm+1
            z2=z2+zlm(lm)*ylmgp(lm,ig)
          ENDDO 
          t3=t3*t2
          z1=z1+t3*z2
        ENDDO 
        z2=jlgprmt(lnpsd,ig,is)*conjg(sfacgp(ig,ias))
        zvclir(jg)=zvclir(jg)+z1*z2
      ELSE 
        t2=y00/factnm(2*lnpsd+1,2)
        zvclir(jg)=zvclir(jg)+t2*zlm(1)
      ENDIF 
    ENDDO 
  ENDDO 
  
  ! solve Poisson's equation in G+p-space for the pseudocharge
  DO ig=1,ngp
    jg=igf(ig)
    zvclir(jg)=gclgp(ig)*zvclir(jg)
  ENDDO 
  
  ! match potentials at muffin-tin boundary by adding homogeneous solution
  DO ias=1,natmtot
    is=idxis(ias)
  ! find the spherical harmonic expansion of the interstitial potential at the
  ! muffin-tin radius
    zlm(:)=0.d0
    DO ig=1,ngp
      jg=igf(ig)
      z1=fourpi*zvclir(jg)*sfacgp(ig,ias)
      lm=0
      DO l=0,lmaxo
        z2=jlgprmt(l,ig,is)*z1*zil(l)
        DO m=-l,l
          lm=lm+1
          zlm(lm)=zlm(lm)+z2*conjg(ylmgp(lm,ig))
        ENDDO 
      ENDDO 
    ENDDO 
  ! calculate the homogenous solution
    i=np(is)-lmmaxo
    lm=0
    DO l=0,lmaxi
      t1=1.d0/rmtl(l,is)
      DO m=-l,l
        lm=lm+1
        i=i+1
        z1=t1*(zlm(lm)-zvclmt(i,ias))
        j=lm
        DO ir=1,nri(is)
          zhmt(j)=z1*rl(ir,l,is)
          j=j+lmmaxi
        ENDDO 
        DO ir=nri(is)+1,nr(is)
          zhmt(j)=z1*rl(ir,l,is)
          j=j+lmmaxo
        ENDDO 
      ENDDO 
    ENDDO 
    DO l=lmaxi+1,lmaxo
      t1=1.d0/rmtl(l,is)
      DO m=-l,l
        lm=lm+1
        i=i+1
        z1=t1*(zlm(lm)-zvclmt(i,ias))
        j=npi(is)+lm
        DO ir=nri(is)+1,nr(is)
          zhmt(j)=z1*rl(ir,l,is)
          j=j+lmmaxo
        ENDDO 
      ENDDO 
    ENDDO 
    zvclmt(1:np(is),ias)=zvclmt(1:np(is),ias)+zhmt(1:np(is))
  
    ! store the nuclear potential without the self-term for the phonon dynamical
    ! matrix calculation
    !IF(tphdyn) THEN 
    !  IF(ias.eq.iasph) zvnmt(1:np(is))=zhmt(1:np(is))
    !ENDIF 
  
  ENDDO 

  ! Fourier transform interstitial potential to real-space
  CALL zfftifc(3,ngdg,1,zvclir)
  RETURN 
END SUBROUTINE 
