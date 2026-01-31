SUBROUTINE my_energy()
  USE modmain
  USE moddftu
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
  REAL(8), external :: rfinp
  COMPLEX(8), external :: zdotc
  
  ! coupling constant of the external field (g_e/4c)
  cb = gfacte/(4.d0*solsc)

  !-----------------------------------------------!
  !     exchange-correlation potential energy     !
  !-----------------------------------------------!
  engyvxc = rfinp(rhomt, rhoir, vxcmt, vxcir)
  write(*,*) 'engyvxc = ', engyvxc

  !-----------------------------------------------------!
  !     exchange-correlation effective field energy     !
  !-----------------------------------------------------!
  engybxc = 0.d0
  DO idm = 1,ndmag
    engybxc = engybxc + rfinp(magmt(:,:,idm), magir(:,idm), bxcmt(:,:,idm), bxcir(:,idm))
  ENDDO 

  !------------------------------------------!
  !     external magnetic field energies     !
  !------------------------------------------!
  engybext = 0.d0
  DO idm = 1,ndmag
    IF(ncmag) THEN 
      jdm = idm
    else
      jdm = 3
    ENDIF 
    ! energy of physical global field
    engybext = engybext + cb*momtot(idm)*bfieldc(jdm)
  ENDDO 

  !----------------------------------!
  !     Coulomb potential energy     !
  !----------------------------------!
  engyvcl = rfinp(rhomt, rhoir, vclmt, vclir)

  !-----------------------!
  !     Madelung term     !
  !-----------------------!
  engymad = 0.d0
  DO ias = 1,natmtot
    is = idxis(ias)
    engymad = engymad + 0.5d0*spzn(is)*(vclmt(1,ias) - vcln(1,is))*y00
  ENDDO 

  !---------------------------------------------!
  !     electron-nuclear interaction energy     !
  !---------------------------------------------!
  engyen = 2.d0*(engymad - engynn)

  !------------------------!
  !     Hartree energy     !
  !------------------------!
  engyhar = 0.5d0*(engyvcl - engyen)

  !------------------------!
  !     Coulomb energy     !
  !------------------------!
  engycl = engynn + engyen + engyhar
  
  !-------------------------!
  !     exchange energy     !
  !-------------------------!
  IF((xctype(1) < 0) .or. (task == 5)) THEN 
    ! exact exchange for OEP-EXX or Hartree-Fock on last self-consistent loop
    IF(tlast) THEN 
      CALL exxengy()
      ! mix exact and DFT exchange energies for hybrid functionals
      IF(hybrid) THEN 
        engyx = engyx*hybridc
        engyx = engyx + rfinp(rhomt,rhoir,exmt,exir)
      ENDIF 
    ELSE 
      engyx = 0.d0
    ENDIF 
  ELSE
    ! exchange energy from the density
    engyx = rfinp(rhomt, rhoir, exmt, exir)
  ENDIF 

  !----------------------------!
  !     correlation energy     !
  !----------------------------!
  IF(task==5) THEN 
    IF(hybrid) THEN 
      ! fraction of DFT correlation energy for hybrid functionals
      engyc = rfinp(rhomt,rhoir,ecmt,ecir)
    ELSE 
      ! zero correlation energy for pure Hartree-Fock
      engyc = 0.d0
    ENDIF 
  ELSE 
    ! correlation energy from the density
    engyc = rfinp(rhomt, rhoir, ecmt, ecir)
  ENDIF 

  !----------------------!
  !     DFT+U energy     !
  !----------------------!
  engydu =0.d0
  if (dftu /= 0) then
    do i = 1,ndftu
      is = idftu(1,i)
      do ia = 1,natoms(is)
        engydu = engydu + engyadu(ia,i)
      enddo
    enddo
  endif


  !----------------------------!
  !     sum of eigenvalues     !
  !----------------------------!
  ! core eigenvalues
  evalsum = 0.d0
  DO ias = 1,natmtot
    is = idxis(ias)
    DO ist = 1,nstsp(is)
      IF(spcore(ist,is)) then
        evalsum = evalsum + occcr(ist,ias)*evalcr(ist,ias)
      ENDIF
    ENDDO 
  ENDDO 
  ! valence eigenvalues
  DO ik = 1,nkpt
    DO ist = 1,nstsv
      evalsum = evalsum + wkpt(ik)*occsv(ist,ik)*evalsv(ist,ik)
    ENDDO 
  ENDDO 

  !------------------------!
  !     kinetic energy     !
  !------------------------!
  ! core electron kinetic energy
  CALL energykncr()

  ! total electron kinetic energy
  IF(task == 5) THEN 
    ! Hartree-Fock case
    engykn = engykncr
    ! kinetic energy from valence states
    ALLOCATE(evecsv(nstsv,nstsv),kmat(nstsv,nstsv),c(nstsv,nstsv))
    DO ik = 1,nkpt
      CALL getevecsv(filext,ik,vkl(:,ik),evecsv)
      CALL getkmat(ik,kmat)
      CALL zgemm('N','N',nstsv,nstsv,nstsv,zone,kmat,nstsv,evecsv,nstsv,zzero,c,nstsv)
      DO ist=1,nstsv
        z1 = zdotc(nstsv,evecsv(:,ist),1,c(:,ist),1)
        engykn = engykn + wkpt(ik)*occsv(ist,ik)*dble(z1)
      ENDDO 
    ENDDO 
    DEALLOCATE(evecsv,kmat,c)
  ELSE 
    ! Kohn-Sham case
    ALLOCATE(rfmt(npmtmax,natmtot))
    ! remove magnetic field contribution
    sum = 0.d0
    DO idm = 1,ndmag
      DO ias = 1,natmtot
        is = idxis(ias)
        CALL rfsht(nrcmt(is),nrcmti(is),bsmt(:,ias,idm),rfmt(:,ias))
      ENDDO 
      CALL rfmtctof(rfmt)
      sum = sum + rfinp(magmt(:,:,idm),magir(:,idm),rfmt,bsir(:,idm))
    ENDDO 
    write(*,*) 'ss_mag_field = ', sum
    !
    ! remove integral of w_xc times tau for generalised Kohn-Sham meta-GGA
    IF(xcgrad == 4) THEN 
      DO ispn=1,nspinor
        DO ias=1,natmtot
          is=idxis(ias)
          np=npmt(is)
          rfmt(1:np,ias) = taumt(1:np,ias,ispn) - taucr(1:np,ias,ispn)
        ENDDO 
        sum = sum + rfinp(rfmt,tauir,wxcmt,wxcir)
      ENDDO 
    ENDIF 
    ! remove fixed tensor moment potential matrix contribution
    if (ftmtype /= 0) then
      n2=(lmmaxdm*nspinor)**2
      do ias=1,natmtot
        z1=zdotc(n2,dmatmt(:,:,:,:,ias),1,vmftm(:,:,:,:,ias),1)
        sum=sum+dble(z1)
      end do
    end if
    engykn = evalsum - engyvcl - engyvxc - sum
    DEALLOCATE(rfmt)
  ENDIF 
  
  !-------------------------------!
  !     entropic contribution     !
  !-------------------------------!
  entrpy = 0.d0
  engyts = 0.d0
  ! non-zero only for the Fermi-Dirac smearing function
  IF(stype == 3) THEN 
    sum = 0.d0
    DO ik=1,nkpt
      DO ist=1,nstsv
        f=occsv(ist,ik)/occmax
        IF((f > 0.d0) .and. (f < 1.d0)) THEN 
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
  engytot = engykn + 0.5d0*engyvcl + engymad + engyx + engyc + engyts
  
  ! add the DFT+U correction if required
  if (dftu /= 0) engytot = engytot + engydu

  ! write total energy
  WRITE(*,*) 'total energy', engytot
  RETURN 
END SUBROUTINE 
