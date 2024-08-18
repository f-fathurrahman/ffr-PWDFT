! evecsv argument is removed

! input:
! - ik: kpoint index
! - evalfv: eigenvalues
! - evecfv: eigenvectors
SUBROUTINE my_eveqn( ik, evalfv, evecfv )
  USE m_atoms, ONLY: natmtot
  USE m_gkvectors, ONLY: ngk, ngkmax, igkig, vgkc, gkc, sfacgk
  USE m_hamiltonian, ONLY: nmatmax, nmat
  USE m_states, ONLY: nstfv
  USE m_kpoints, ONLY: vkc
  USE m_muffin_tins, ONLY: lmmaxapw
  USE m_apwlo, ONLY: apwordmax
  USE m_spin, ONLY: nspnfv
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ik
  REAL(8), intent(out) :: evalfv(nstfv,nspnfv)
  COMPLEX(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv)
  ! local variables
  INTEGER jspn
  ! allocatable arrays
  COMPLEX(8), ALLOCATABLE :: apwalm(:,:,:,:,:)
  ALLOCATE(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  
  write(*,*)
  write(*,*) '**** ENTER my_eveqn ****'
  write(*,*)

  ! loop over first-variational spins (nspnfv=2 for spin-spirals only)
  DO jspn = 1,nspnfv
    write(*,'(1x,A,I5,A,I5)') 'ik = ', ik, ' jspn = ', jspn
    ! find the matching coefficients
    CALL match( ngk(jspn,ik), vgkc(:,:,jspn,ik), gkc(:,jspn,ik), &
                sfacgk(:,:,jspn,ik), apwalm(:,:,:,:,jspn) )
    !
    ! the important variable here is apwalm, which will be passed to eveqnfv
    !
    !
    ! solve the first-variational eigenvalue equation
    !
    ! XXX: tefvit is disabled
    ! directly
    CALL my_eveqnfv( nmat(jspn,ik), ngk(jspn,ik), igkig(:,jspn,ik), vkc(:,ik), &
                  vgkc(:,:,jspn,ik), apwalm(:,:,:,:,jspn), evalfv(:,jspn), evecfv(:,:,jspn))
    ! NOTE: ik index is not passed directly
  ENDDO 
  
  ! XXX: spinsprl and second variational stuffs are removed

  DEALLOCATE(apwalm)

  write(*,*)
  write(*,*) '**** EXIT my_eveqn ****'
  write(*,*)


  RETURN
END SUBROUTINE 
