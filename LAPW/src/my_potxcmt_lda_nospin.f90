! a simplified version of potxcmt for LDA and nospin case
SUBROUTINE my_potxcmt_lda_nospin( &
  tsh, ias, xctype_, &
  rhomt_, &
  exmt_, ecmt_, vxcmt_)
  !
  USE m_atoms, ONLY: natmtot, idxis
  USE m_spin, ONLY: spinpol
  USE m_muffin_tins, ONLY: npmtmax, npmt, nrmt, nrmti
  USE m_density_pot_xc, ONLY: xcgrad
  USE m_states, ONLY: swidth
  USE m_oep_hf, ONLY: hybridc, hybrid
  USE modxcifc, ONLY: xcifc
  !
  IMPLICIT NONE 
  ! arguments
  logical, intent(in) :: tsh
  INTEGER, intent(in) :: ias,xctype_(3)
  REAL(8), intent(in) :: rhomt_(npmtmax,natmtot)
  REAL(8), intent(out) :: exmt_(npmtmax,natmtot), ecmt_(npmtmax,natmtot)
  REAL(8), intent(out) :: vxcmt_(npmtmax,natmtot)
  ! local variables
  INTEGER is
  INTEGER nr,nri,n
  REAL(8) t1
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rho(:)
  REAL(8), ALLOCATABLE :: ex(:),ec(:),vxc(:)
  REAL(8), ALLOCATABLE :: vx(:)
  REAL(8), ALLOCATABLE :: vc(:)

  is = idxis(ias)
  n = npmt(is)

  ! allocate local arrays
  ALLOCATE(rho(n),ex(n),ec(n),vxc(n))

  if( xcgrad /= 0 ) then
    stop 'my_potxcmt_lda_nospin: xcgrad /= 0 is not supported in this subroutine'
  endif

  IF(spinpol) THEN 
    stop 'my_potxcmt_lda_nospin: spinpol == true is not supported in this subroutine'
  ENDIF 

  ! No spinpol
  ALLOCATE(vx(n), vc(n))

  nr = nrmt(is)
  nri = nrmti(is)

  ! tsh seems to be fixed to .true.
  IF(tsh) THEN 
    ! convert the density to spherical coordinates
    CALL my_rbsht(nr, nri, rhomt_(:,ias), rho)
  ELSE 
    rho(1:n) = rhomt_(1:n,ias)
  ENDIF 

  ! convert tau to spherical coordinates if required
  ! ... Not applicable here

  ! LDA, spin-unpolarized case only
  CALL xcifc(xctype_, n=n, tempa=swidth, rho=rho, ex=ex, ec=ec, vx=vx, vc=vc)

  ! hybrid functionals
  IF(hybrid) THEN 
    t1 = 1.d0 - hybridc
    ! scale exchange part of energy
    ex(:) = t1*ex(:)
    ! scale exchange part of potential
    vxc(:) = t1*vx(:) + vc(:)
  ELSE 
    vxc(:) = vx(:) + vc(:)
  ENDIF 


  write(*,*)
  write(*,*) 'my_potxcmt_lda_nospin: tsh = ', tsh 
  IF(tsh) THEN 
    ! convert exchange and correlation energy densities to spherical harmonics
    CALL my_rfsht(nr, nri, ex, exmt_(:,ias))
    CALL my_rfsht(nr, nri, ec, ecmt_(:,ias))
    ! convert exchange-correlation potential to spherical harmonics
    CALL my_rfsht(nr, nri, vxc, vxcmt_(:,ias))
    write(*,*) 'my_potxcmt: shape(exmt)  = ', shape(exmt_)
    write(*,*) 'my_potxcmt: shape(ecmt)  = ', shape(ecmt_)
    write(*,*) 'my_potxcmt: shape(vxcmt) = ', shape(vxcmt_)
  ELSE 
    exmt_(1:n,ias) = ex(1:n)
    ecmt_(1:n,ias) = ec(1:n)
    vxcmt_(1:n,ias) = vxc(1:n)
  ENDIF 

  ! Deallocate memory
  DEALLOCATE(rho, ex, ec, vxc)
  DEALLOCATE(vx,vc)

  RETURN 
END SUBROUTINE 

