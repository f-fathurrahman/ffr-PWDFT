! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point number (in,integer)
!   evalfv : first-variational eigenvalues (out,real(nstfv))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (out,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Solves the first- and second-variational eigenvalue equations. See routines
!   {\tt match}, {\tt eveqnfv}, {\tt eveqnss} and {\tt eveqnsv}.
SUBROUTINE eveqn(ik,evalfv,evecfv,evecsv)
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
  COMPLEX(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
  ! local variables
  INTEGER jspn
  ! allocatable arrays
  COMPLEX(8), ALLOCATABLE :: apwalm(:,:,:,:,:)
  ALLOCATE(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  
  ! loop over first-variational spins (nspnfv=2 for spin-spirals only)
  DO jspn=1,nspnfv
    ! find the matching coefficients
    CALL match(ngk(jspn,ik),vgkc(:,:,jspn,ik),gkc(:,jspn,ik), &
     sfacgk(:,:,jspn,ik),apwalm(:,:,:,:,jspn))
    ! solve the first-variational eigenvalue equation
    IF(tefvit) THEN 
      ! iteratively
      CALL eveqnit(nmat(jspn,ik),ngk(jspn,ik),igkig(:,jspn,ik),vkl(:,ik), &
       vgkl(:,:,jspn,ik),vgkc(:,:,jspn,ik),apwalm(:,:,:,:,jspn),evalfv(:,jspn), &
       evecfv(:,:,jspn))
    ELSE 
      ! directly
      CALL eveqnfv(nmat(jspn,ik),ngk(jspn,ik),igkig(:,jspn,ik),vkc(:,ik), &
       vgkc(:,:,jspn,ik),apwalm(:,:,:,:,jspn),evalfv(:,jspn),evecfv(:,:,jspn))
    ENDIF 
  ENDDO 
  
  IF(spinsprl) THEN 
    ! solve the spin-spiral second-variational eigenvalue equation
    CALL eveqnss(ngk(:,ik),igkig(:,:,ik),apwalm,evalfv,evecfv,evalsv(:,ik),evecsv)
  ELSE 
    ! solve the second-variational eigenvalue equation
    CALL eveqnsv(ngk(1,ik),igkig(:,1,ik),vgkc(:,:,1,ik),apwalm,evalfv,evecfv, &
     evalsv(:,ik),evecsv)
  ENDIF 
  DEALLOCATE(apwalm)
  RETURN
END SUBROUTINE 
