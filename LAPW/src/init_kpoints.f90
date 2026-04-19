SUBROUTINE init_kpoints()
  USE m_constants, ONLY: twopi
  USE m_symmetry, ONLY: lsplsymc, lspnsymc, tv0symc, nsymlat, nsymcrys, symlat
  USE m_kpoints, ONLY: ngridk, kptboxl, wkpt, wkptnr, vkl, vkloff, vkc, symkpt, reducek, &
                 radkpt, nsymkpt, nkptnr, nkpt, ivk, ivkik, ivkiknr, autokpt, ndspem, vklem, deltaem
  USE m_misc, ONLY: task
  USE m_lattice, ONLY: epslat, bvec, binv
  USE m_atoms, ONLY: molecule
  USE m_spin, ONLY: spinpol
  USE m_plotting
  IMPLICIT NONE 
  LOGICAL :: lsym(48)
  INTEGER :: isym, ik, i1, i2, i3
  REAL(8) :: t1, vc(3), vl(3)

  !---------------------!
  !     k-point set     !
  !---------------------!
  !
  ! check if the system is an isolated molecule
  IF( molecule ) THEN 
    ngridk(:) = 1
    vkloff(:) = 0.d0
    autokpt = .false.
  ENDIF 
  !
  ! store the point group symmetries for reducing the k-point set
  IF( reducek == 0 ) THEN 
    nsymkpt = 1
    symkpt(:,:,1) = symlat(:,:,1)
  ELSE 
    lsym(:) = .false.
    DO isym = 1,nsymcrys
      IF( reducek == 2 ) THEN 
        ! check symmetry is symmorphic
        IF( .not. tv0symc(isym)) GOTO 10
        ! check also that the spin rotation is the same as the spatial rotation
        IF( spinpol ) THEN 
          IF( lspnsymc(isym) /= lsplsymc(isym) ) GOTO 10
        ENDIF 
      ENDIF 
      lsym(lsplsymc(isym)) = .true.
      10 CONTINUE 
    ENDDO 
    nsymkpt = 0
    DO isym = 1,nsymlat
      IF( lsym(isym) ) THEN 
        nsymkpt = nsymkpt + 1
        symkpt(:,:,nsymkpt)=symlat(:,:,isym)
      ENDIF 
    ENDDO 
  ENDIF 

  if (any(task == [20,21,22,23])) then
    ! generate k-points along a path for band structure plots
    call plotpt1d(bvec, nvp1d, npp1d, vvlp1d,vplp1d,dvp1d,dpp1d)
    nkpt = npp1d
    !
    if(allocated(vkl)) deallocate(vkl)
    allocate(vkl(3,nkpt))
    !
    if (allocated(vkc)) deallocate(vkc)
    allocate(vkc(3,nkpt))
    do ik=1,nkpt
      vkl(:,ik) = vplp1d(:,ik)
      call r3mv(bvec, vkl(:,ik), vkc(:,ik))
    enddo
    nkptnr = nkpt
  ELSEIF( task == 25 ) THEN 
    ! effective mass calculation
    nkpt = (2*ndspem+1)**3
    !
    if (allocated(ivk)) deallocate(ivk)
    allocate(ivk(3,nkpt))
    !
    if (allocated(vkl)) deallocate(vkl)
    allocate(vkl(3,nkpt))
    !
    if (allocated(vkc)) deallocate(vkc)
    allocate(vkc(3,nkpt))
    ! map vector to [0,1)
    call r3frac(epslat,vklem)
    ik = 0
    do i3 = -ndspem,ndspem
      do i2 = -ndspem,ndspem
        do i1 = -ndspem,ndspem
          ik = ik+1
          ivk(1,ik) = i1; ivk(2,ik)=i2; ivk(3,ik)=i3
          vc(1) = dble(i1); vc(2)=dble(i2); vc(3)=dble(i3)
          vc(:) = vc(:)*deltaem
          call r3mv(binv,vc,vl)
          vkl(:,ik) = vklem(:) + vl(:)
          call r3mv(bvec,vkl(:,ik),vkc(:,ik))
        enddo
      enddo
    enddo
    nkptnr=nkpt
  ELSE
    ! determine the k-point grid automatically from radkpt if required
    IF( autokpt ) THEN 
      t1 = radkpt/twopi
      ngridk(:) = int(t1*sqrt(bvec(1,:)**2 + bvec(2,:)**2 + bvec(3,:)**2)) + 1
    ENDIF 
    ! set up the default k-point box
    kptboxl(:,0) = vkloff(:)/dble(ngridk(:))
    IF(task == 102) kptboxl(:,0)=0.d0
    !
    kptboxl(:,1) = kptboxl(:,0)
    kptboxl(:,2) = kptboxl(:,0)
    kptboxl(:,3) = kptboxl(:,0)
    kptboxl(1,1) = kptboxl(1,1)+1.d0
    kptboxl(2,2) = kptboxl(2,2)+1.d0
    kptboxl(3,3) = kptboxl(3,3)+1.d0
    ! k-point set and box for Fermi surface plots
    IF( any(task==[100,101,102] )) THEN 
      ngridk(:) = np3d(:)
      IF( task /= 102 ) kptboxl(:,:) = vclp3d(:,:)
    ENDIF 
    !
    ! allocate the k-point set arrays
    IF( allocated(ivkik) ) DEALLOCATE(ivkik)
    ALLOCATE( ivkik(0:ngridk(1)-1, 0:ngridk(2)-1, 0:ngridk(3)-1) )
    !
    IF( allocated(ivkiknr) ) DEALLOCATE(ivkiknr)
    ALLOCATE( ivkiknr(0:ngridk(1)-1, 0:ngridk(2)-1, 0:ngridk(3)-1) )
    !
    nkptnr = ngridk(1)*ngridk(2)*ngridk(3)
    !
    IF( allocated(ivk)) DEALLOCATE(ivk)
    ALLOCATE( ivk(3,nkptnr) )
    IF( allocated(vkl) ) DEALLOCATE(vkl)
    ALLOCATE( vkl(3,nkptnr) )
    IF( allocated(vkc) ) DEALLOCATE(vkc)
    ALLOCATE( vkc(3,nkptnr))
    IF( allocated(wkpt) ) DEALLOCATE(wkpt)
    ALLOCATE( wkpt(nkptnr) )
    ! generate the k-point set
    CALL genppts( .false., nsymkpt, symkpt, ngridk, nkptnr, &
                  epslat,bvec,kptboxl,nkpt, &
                  ivkik, ivkiknr, ivk, vkl, vkc, wkpt, wkptnr )
  ENDIF 

END SUBROUTINE 
