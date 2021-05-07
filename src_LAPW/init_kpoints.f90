SUBROUTINE init_kpoints()
  
  USE modmain
  IMPLICIT NONE 
  LOGICAL :: lsym(48)
  INTEGER :: isym
  REAL(8) :: t1

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

  if (any(task.eq.[20,21,22,23])) then
    ! generate k-points along a path for band structure plots
    ! XXX Disabled
    WRITE(*,*) 'task = ', task
    WRITE(*,*) 'disabled ...'
    STOP 
  ELSEIF( task == 25 ) THEN 
    ! effective mass calculation
    ! XXX Disabled
    WRITE(*,*) 'task = ', task
    WRITE(*,*) 'disabled ...'
    STOP 
  ELSE
    ! determine the k-point grid automatically from radkpt if required
    IF( autokpt ) THEN 
      t1 = radkpt/twopi
      ngridk(:) = int(t1*sqrt(bvec(1,:)**2 + bvec(2,:)**2 + bvec(3,:)**2)) + 1
    ENDIF 
    ! set up the default k-point box
    kptboxl(:,0) = vkloff(:)/dble(ngridk(:))
    IF(task.eq.102) kptboxl(:,0)=0.d0
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
