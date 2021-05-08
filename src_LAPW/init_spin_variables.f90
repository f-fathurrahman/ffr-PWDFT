!-------------------------------
SUBROUTINE init_spin_variables()
!------------------------------
  USE modmain, ONLY: &
           spinsprl, ssdph, natoms, natmtot, xctype, jspnfv, xcgrad, xcdescr, xcspin, &
           tevecsv, task, tefvit, spinpol, spincore, spinorb, spcpl, occmax, reducebf, &
           nspinor, nspnfv, nspecies, ncmag, ndmag, nosource, mixtype, &
           hybrid, hybridc, epslat, cmagz, fsmtype, bfcmt0, bfcmt, bfieldc0, bfieldc, &
           bfsmc, bfsmcmt
  USE modxcifc, ONLY: getxcdata

  IMPLICIT NONE 
  INTEGER :: i, ia, is

  !------------------------!
  !     spin variables     !
  !------------------------!
  IF( spinsprl ) THEN 
    spinpol = .true.
    spinorb = .false.
    !
    IF( any(task == [51,52,53,61,62,63,700,701]) ) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(init0): spin-spirals do not work with task ",I4)') task
      WRITE(*,*)
      STOP 
    ENDIF 
    !
    IF( xctype(1) .lt. 0 ) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(init0): spin-spirals do not work with the OEP method")')
      WRITE(*,*)
      STOP 
    ENDIF 
    !
  ENDIF 
  
  ! de-phasing required only for spin-spirals
  IF( .not. spinsprl ) ssdph = .false.
  
  ! spin-orbit coupling, fixed spin moment, spin spirals or spin-polarised cores
  ! requires a spin-polarised calculation
  IF( (spinorb) .or. (fsmtype /= 0) .or. (spinsprl) .or. (spincore) ) spinpol = .true.
  
  ! number of spinor components and maximum allowed occupancy
  IF( spinpol ) THEN 
    nspinor = 2
    occmax = 1.d0
  ELSE 
    nspinor = 1
    occmax = 2.d0
  ENDIF 
  
  ! number of spin-dependent first-variational functions per state and map from
  ! second- to first-variational spin index
  IF( spinsprl ) THEN 
    nspnfv = 2
    jspnfv(1) = 1
    jspnfv(2) = 2
  ELSE 
    nspnfv = 1
    jspnfv(1) = 1
    jspnfv(2) = 1
  ENDIF 
  
  ! no calculation of second-variational eigenvectors by default
  tevecsv = .false.
  
  ! spin-polarised calculations require second-variational eigenvectors
  IF( spinpol ) tevecsv = .true.
  
  ! Hartree-Fock/RDMFT/TDDFT/GW requires second-variational eigenvectors
  IF( any( task == [5,10,170,300,460,461,600,620,630]) ) THEN 
    tevecsv = .true.
  ENDIF 
  
  ! get exchange-correlation functional data
  CALL getxcdata(xctype, xcdescr, xcspin, xcgrad, hybrid, hybridc)
  IF( (spinpol) .and. (xcspin == 0) ) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(init0): requested spin-polarised run with spin-unpolarised")')
    write(*,'(" exchange-correlation functional")')
    write(*,*)
    STOP 
  ENDIF 
  
  ! check for collinearity in the z-direction and set the dimension of the
  ! magnetisation and exchange-correlation vector fields
  IF( spinpol) THEN 
    ndmag = 1
    IF( (abs(bfieldc0(1)) > epslat) .or. (abs(bfieldc0(2)) > epslat) ) ndmag=3
    DO is = 1,nspecies
      DO ia = 1,natoms(is)
        IF( (abs(bfcmt0(1,ia,is)) > epslat) .or. &
            (abs(bfcmt0(2,ia,is)) > epslat)) ndmag=3
      ENDDO 
    ENDDO 
    ! spin-orbit coupling is non-collinear in general
    IF( spinorb ) ndmag=3
    ! source-free fields and spin-spirals must be non-collinear
    IF( (nosource) .or. (spinsprl) ) THEN 
      ndmag = 3
      cmagz = .false.
    ENDIF 
    ! force collinear magnetism along the z-axis if required
    IF( cmagz ) ndmag=1
  ELSE 
    ndmag = 0
  ENDIF 
  
  ! set the non-collinear flag
  IF( ndmag == 3 ) THEN 
    ncmag = .true.
  ELSE 
    ncmag = .false.
  ENDIF 

  ! check for meta-GGA with non-collinearity
  IF( ( (xcgrad == 3) .or. (xcgrad == 4)) .and. ncmag) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(init0): meta-GGA is not valid for non-collinear magnetism")')
    WRITE(*,*)
    STOP
  ENDIF 

  IF( ncmag .or. spinorb ) THEN 
    ! spins are coupled in the second-variational Hamiltonian
    spcpl = .true.
  ELSE 
    ! spins are decoupled
    spcpl = .false.
  ENDIF 

  ! spin-polarised cores
  IF( .not. spinpol) spincore = .false.
  IF(fsmtype /= 0) THEN 
    ! set fixed spin moment effective field to zero
    bfsmc(:) = 0.d0
    ! set muffin-tin FSM fields to zero
    IF( allocated(bfsmcmt) ) DEALLOCATE(bfsmcmt)
    ALLOCATE( bfsmcmt(3,natmtot) )
    bfsmcmt(:,:) = 0.d0
    IF( mixtype /= 1 ) THEN 
      WRITE(*,*)
      WRITE(*,'("Info(init0): mixtype changed to 1 for FSM calculation")')
    ENDIF 
    mixtype = 1
  ENDIF 

  ! XXX: number of independent spin components of the f_xc spin tensor
  ! XXX: nxfxc is not set (for TDDFT)

  ! set the magnetic fields to the initial values
  bfieldc(:) = bfieldc0(:)
  bfcmt(:,:,:) = bfcmt0(:,:,:)
  
  ! if reducebf < 1 then reduce the external magnetic fields immediately for
  ! non-self-consistent calculations or resumptions
  IF( reducebf < 1.d0-1.d-4 ) THEN 
    IF( all( task /= [0,1,2,3,5,28,200,201,350,351,360,630]) ) THEN 
      bfieldc(:) = 0.d0
      bfcmt(:,:,:) = 0.d0
    ENDIF 
  ENDIF 

  ! set the fixed tensor moment spatial and spin rotation matrices equal for the
  ! case of spin-orbit coupling; parity for spin is ignored by rotdmat
  IF( spinorb ) THEN 
    WRITE(*,*) 'Not supported yet'
    STOP
  ENDIF 
  
  ! write to VARIABLES.OUT
  WRITE(*,*) 'nspinor = ', nspinor
  WRITE(*,*) 'ndmag   = ', ndmag

END SUBROUTINE 
