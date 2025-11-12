! evecsv argument is removed

! input:
! - ik: kpoint index
! - evalfv: eigenvalues
! - evecfv: eigenvectors
SUBROUTINE my_eveqn( ik, evalfv, evecfv, evecsv )
  USE m_atoms, ONLY: natmtot
  USE m_gkvectors, ONLY: ngk, ngkmax, igkig, vgkc, vgkl, gkc, sfacgk
  USE m_hamiltonian, ONLY: tefvit, nmatmax, nmat
  USE m_states, ONLY: nstfv, nstsv, evalsv
  USE m_kpoints, ONLY: vkc, vkl
  USE m_muffin_tins, ONLY: lmmaxapw
  USE m_apwlo, ONLY: apwordmax
  USE m_spin, ONLY: spinsprl, nspnfv
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ik
  REAL(8), intent(out) :: evalfv(nstfv,nspnfv)
  COMPLEX(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv), evecsv(nstsv,nstsv)
  ! local variables
  INTEGER jspn
  ! allocatable arrays
  COMPLEX(8), ALLOCATABLE :: apwalm(:,:,:,:,:)
  ALLOCATE(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  
  write(*,*)
  write(*,*) '<div> ENTER my_eveqn: ik = ', ik
  write(*,*)

  ! loop over first-variational spins (nspnfv=2 for spin-spirals only)
  ! for spinpol nspnfv = 1 (why?)
  DO jspn=1,nspnfv
    !
    write(*,'(1x,A,I5,A,I5)') 'Loop jspn = ', jspn, ' until nspnfv = ', nspnfv
    !
    ! find the matching coefficients
    CALL match(ngk(jspn,ik),vgkc(:,:,jspn,ik),gkc(:,jspn,ik), &
               sfacgk(:,:,jspn,ik),apwalm(:,:,:,:,jspn))
    !
    ! solve the first-variational eigenvalue equation
    !
    IF( tefvit ) THEN 
      ! iteratively
      CALL eveqnit( nmat(jspn,ik), ngk(jspn,ik), igkig(:,jspn,ik), vkl(:,ik), &
                    vgkl(:,:,jspn,ik), vgkc(:,:,jspn,ik), apwalm(:,:,:,:,jspn),evalfv(:,jspn), &
                    evecfv(:,:,jspn))
    ELSE 
      ! directly (for first-variational)
      CALL my_eveqnfv( nmat(jspn,ik), ngk(jspn,ik), igkig(:,jspn,ik), vkc(:,ik), &
                       vgkc(:,:,jspn,ik), apwalm(:,:,:,:,jspn), &
                       evalfv(:,jspn),evecfv(:,:,jspn) )
      ! NOTE: ik index is not passed directly
    ENDIF 
  ENDDO   

  ! Now solve second-variational
  IF(spinsprl) THEN 
    ! solve the spin-spiral second-variational eigenvalue equation
    CALL eveqnss(ngk(:,ik),igkig(:,:,ik),apwalm,evalfv,evecfv,evalsv(:,ik),evecsv)
  ELSE 
    ! solve the second-variational eigenvalue equation
    CALL my_eveqnsv( ngk(1,ik), igkig(:,1,ik), vgkc(:,:,1,ik), apwalm, &
                     evalfv, evecfv, &
                     evalsv(:,ik), evecsv)
  ENDIF 

  write(*,*)
  write(*,*) '</div> EXIT my_eveqn'
  write(*,*)

  DEALLOCATE(apwalm)

  RETURN
END SUBROUTINE 
