SUBROUTINE energy
  USE modmain
! !DESCRIPTION:
!   Computes the total energy and its individual contributions. The kinetic
!   energy is given by
!   $$ T_s=\sum_i n_i\epsilon_i-\int\rho({\bf r})[v_{\rm C}({\bf r})
!    +v_{\rm xc}({\bf r})]d{\bf r}-\int {\bf m}({\bf r})\cdot
!    ({\bf B}_{\rm xc}({\bf r})+{\bf B}_{\rm ext}({\bf r}))d{\bf r}, $$
!   where $n_i$ are the occupancies and $\epsilon_i$ are the eigenvalues of both
!   the core and valence states; $\rho$ is the density; ${\bf m}$ is the
!   magnetisation density; $v_{\rm C}$ is the Coulomb potential; $v_{\rm xc}$
!   and ${\bf B}_{\rm xc}$ are the exchange-correlation potential and magnetic
!   field, respectively; and ${\bf B}_{\rm ext}$ is the external magnetic field.
!   The Hartree, electron-nuclear and nuclear-nuclear electrostatic energies are
!   combined into the Coulomb energy:
!   \begin{align*}
!    E_{\rm C}&=E_{\rm H}+E_{\rm en}+E_{\rm nn} \\
!             &=\frac{1}{2}V_{\rm C}+E_{\rm Mad},
!   \end{align*}
!   where
!   $$ V_{\rm C}=\int\rho({\bf r})v_{\rm C}({\bf r})d{\bf r} $$
!   is the Coulomb potential energy. The Madelung energy is given by
!   $$ E_{\rm Mad}=\frac{1}{2}\sum_{\alpha}z_{\alpha}R_{\alpha}, $$
!   where
!   $$ R_{\alpha}=\lim_{r\rightarrow 0}\left(v^{\rm C}_{\alpha;00}(r)Y_{00}
!    +\frac{z_{\alpha}}{r}\right) $$
!   for atom $\alpha$, with $v^{\rm C}_{\alpha;00}$ being the $l=0$ component of
!   the spherical harmonic expansion of $v_{\rm C}$ in the muffin-tin, and
!   $z_{\alpha}$ is the nuclear charge. Using the nuclear-nuclear energy
!   determined at the start of the calculation, the electron-nuclear and Hartree
!   energies can be isolated with
!   $$ E_{\rm en}=2\left(E_{\rm Mad}-E_{\rm nn}\right) $$
!   and
!   $$ E_{\rm H}=\frac{1}{2}(E_{\rm C}-E_{\rm en}). $$
!   Finally, the total energy is
!   $$ E=T_s+E_{\rm C}+E_{\rm xc}, $$
!   where $E_{\rm xc}$ is obtained either by integrating the
!   exchange-correlation energy density, or in the case of exact exchange, the
!   explicit calculation of the Fock exchange integral. The energy from the
!   external magnetic fields in the muffin-tins, {\tt bfcmt}, is always removed
!   from the total since these fields are non-physical: their field lines DO not
!   close. The energy of the physical external field, {\tt bfieldc}, is also not
!   included in the total because this field, like those in the muffin-tins, is
!   used for breaking spin symmetry and taken to be infintesimal. If this field
!   is intended to be finite, THEN  the associated energy, {\tt engybext}, should
!   be added to the total by hand. See {\tt potxc}, {\tt exxengy} and related
!   subroutines.
  IMPLICIT NONE 
  ! local variables
  INTEGER :: ik,ist,ispn,idm,jdm
  INTEGER :: is,ia,ias,np,n2,i
  REAL(8) :: cb,sum,f
  COMPLEX(8) :: z1
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rfmt(:,:)
  COMPLEX(8), ALLOCATABLE :: evecsv(:,:),kmat(:,:),c(:,:)

  ! external functions
  REAL(8) rfinp
  COMPLEX(8) zdotc
  external rfinp,zdotc
  
  ! coupling constant of the external field (g_e/4c)
  cb=gfacte/(4.d0*solsc)

  !-----------------------------------------------!
  !     exchange-correlation potential energy     !
  !-----------------------------------------------!
  engyvxc=rfinp(rhomt,rhoir,vxcmt,vxcir)

  !-----------------------------------------------------!
  !     exchange-correlation effective field energy     !
  !-----------------------------------------------------!
  engybxc=0.d0
  DO idm=1,ndmag
    engybxc=engybxc+rfinp(magmt(:,:,idm),magir(:,idm),bxcmt(:,:,idm),bxcir(:,idm))
  ENDDO 

  !------------------------------------------!
  !     external magnetic field energies     !
  !------------------------------------------!
  engybext=0.d0
  DO idm=1,ndmag
    IF(ncmag) THEN 
      jdm=idm
    else
      jdm=3
    ENDIF 
    ! energy of physical global field
    engybext=engybext+cb*momtot(idm)*bfieldc(jdm)
  ENDDO 

  !----------------------------------!
  !     Coulomb potential energy     !
  !----------------------------------!
  engyvcl=rfinp(rhomt,rhoir,vclmt,vclir)

  !-----------------------!
  !     Madelung term     !
  !-----------------------!
  engymad=0.d0
  DO ias=1,natmtot
    is=idxis(ias)
    engymad=engymad+0.5d0*spzn(is)*(vclmt(1,ias)-vcln(1,is))*y00
  ENDDO 

  !---------------------------------------------!
  !     electron-nuclear interaction energy     !
  !---------------------------------------------!
  engyen=2.d0*(engymad-engynn)

  !------------------------!
  !     Hartree energy     !
  !------------------------!
  engyhar=0.5d0*(engyvcl-engyen)

  !------------------------!
  !     Coulomb energy     !
  !------------------------!
  engycl=engynn+engyen+engyhar
  
  !-------------------------!
  !     exchange energy     !
  !-------------------------!
  IF((xctype(1).lt.0).or.(task.eq.5)) THEN 
    ! exact exchange for OEP-EXX or Hartree-Fock on last self-consistent loop
    IF(tlast) THEN 
      CALL exxengy()
      ! mix exact and DFT exchange energies for hybrid functionals
      IF(hybrid) THEN 
        engyx=engyx*hybridc
        engyx=engyx+rfinp(rhomt,rhoir,exmt,exir)
      ENDIF 
    ELSE 
      engyx=0.d0
    ENDIF 
  ELSE
    ! exchange energy from the density
    engyx=rfinp(rhomt,rhoir,exmt,exir)
  ENDIF 

  !----------------------------!
  !     correlation energy     !
  !----------------------------!
  IF(task==5) THEN 
    IF(hybrid) THEN 
      ! fraction of DFT correlation energy for hybrid functionals
      engyc=rfinp(rhomt,rhoir,ecmt,ecir)
    ELSE 
      ! zero correlation energy for pure Hartree-Fock
      engyc=0.d0
    ENDIF 
  ELSE 
    ! correlation energy from the density
    engyc=rfinp(rhomt,rhoir,ecmt,ecir)
  ENDIF 

  !----------------------------!
  !     sum of eigenvalues     !
  !----------------------------!
  ! core eigenvalues
  evalsum=0.d0
  DO ias=1,natmtot
    is=idxis(ias)
    DO ist=1,nstsp(is)
      IF(spcore(ist,is)) evalsum=evalsum+occcr(ist,ias)*evalcr(ist,ias)
    ENDDO 
  ENDDO 
  ! valence eigenvalues
  DO ik=1,nkpt
    DO ist=1,nstsv
      evalsum=evalsum+wkpt(ik)*occsv(ist,ik)*evalsv(ist,ik)
    ENDDO 
  ENDDO 

  !------------------------!
  !     kinetic energy     !
  !------------------------!
  ! core electron kinetic energy
  CALL energykncr()

  ! total electron kinetic energy
  IF(task.eq.5) THEN 
    ! Hartree-Fock case
    engykn=engykncr
    ! kinetic energy from valence states
    ALLOCATE(evecsv(nstsv,nstsv),kmat(nstsv,nstsv),c(nstsv,nstsv))
    DO ik=1,nkpt
      CALL getevecsv(filext,ik,vkl(:,ik),evecsv)
      CALL getkmat(ik,kmat)
      CALL zgemm('N','N',nstsv,nstsv,nstsv,zone,kmat,nstsv,evecsv,nstsv,zzero,c, &
       nstsv)
      DO ist=1,nstsv
        z1=zdotc(nstsv,evecsv(:,ist),1,c(:,ist),1)
        engykn=engykn+wkpt(ik)*occsv(ist,ik)*dble(z1)
      ENDDO 
    ENDDO 
    DEALLOCATE(evecsv,kmat,c)
  ELSE 
    ! Kohn-Sham case
    ALLOCATE(rfmt(npmtmax,natmtot))
    ! remove magnetic field contribution
    sum=0.d0
    DO idm=1,ndmag
      DO ias=1,natmtot
        is=idxis(ias)
        CALL rfsht(nrcmt(is),nrcmti(is),bsmt(:,ias,idm),rfmt(:,ias))
      ENDDO 
      CALL rfmtctof(rfmt)
      sum=sum+rfinp(magmt(:,:,idm),magir(:,idm),rfmt,bsir(:,idm))
    ENDDO 
    ! remove integral of w_xc times tau for generalised Kohn-Sham meta-GGA
    IF(xcgrad.eq.4) THEN 
      DO ispn=1,nspinor
        DO ias=1,natmtot
          is=idxis(ias)
          np=npmt(is)
          rfmt(1:np,ias)=taumt(1:np,ias,ispn)-taucr(1:np,ias,ispn)
        ENDDO 
        sum=sum+rfinp(rfmt,tauir,wxcmt,wxcir)
      ENDDO 
    ENDIF 
    engykn=evalsum-engyvcl-engyvxc-sum
    DEALLOCATE(rfmt)
  ENDIF 
  
  !-------------------------------!
  !     entropic contribution     !
  !-------------------------------!
  entrpy=0.d0
  engyts=0.d0
  ! non-zero only for the Fermi-Dirac smearing function
  IF(stype.eq.3) THEN 
    sum=0.d0
    DO ik=1,nkpt
      DO ist=1,nstsv
        f=occsv(ist,ik)/occmax
        IF((f.gt.0.d0).and.(f.lt.1.d0)) THEN 
          sum = sum + wkpt(ik)*(f*log(f)+(1.d0-f)*log(1.d0-f))
        ENDIF 
      ENDDO 
    ENDDO 
    ! entropy
    entrpy = -occmax*kboltz*sum
    ! contribution to free energy
    engyts = -swidth*entrpy/kboltz
  ENDIF 

  !----------------------!
  !     total energy     !
  !----------------------!
  engytot=engykn+0.5d0*engyvcl+engymad+engyx+engyc+engyts
  
  ! write total energy
  WRITE(*,*) 'total energy', engytot
  RETURN 
END SUBROUTINE 
