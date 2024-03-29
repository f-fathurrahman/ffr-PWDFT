! !INPUT/OUTPUT PARAMETERS:
!   nmatp  : order of overlap and Hamiltonian matrices (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vpc    : p-vector in Cartesian coordinates (in,real(3))
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evalfv : first-variational eigenvalues (out,real(nstfv))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
! !DESCRIPTION:
!   Solves the eigenvalue equation,
!   $$ (H-\epsilon O)b=0, $$
!   for the all the first-variational states of the input $p$-point.
SUBROUTINE eveqnfv(nmatp,ngp,igpig,vpc,vgpc,apwalm,evalfv,evecfv)
  USE m_atoms, ONLY: natmtot
  USE m_gkvectors, ONLY: ngkmax
  USE m_timing, ONLY: timemat
  USE m_hamiltonian, ONLY: tefvr, nmatmax
  USE m_states, ONLY: nstfv
  USE m_muffin_tins, ONLY: lmmaxapw
  USE m_apwlo, ONLY: apwordmax
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nmatp,ngp,igpig(ngkmax)
  REAL(8), intent(in) :: vpc(3),vgpc(3,ngkmax)
  COMPLEX(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  REAL(8), intent(out) :: evalfv(nstfv)
  COMPLEX(8), intent(out) :: evecfv(nmatmax,nstfv)
  ! local variables
  REAL(8) ts0,ts1
  ! allocatable arrays
  COMPLEX(8), allocatable :: h(:,:),o(:,:)
  
  !-----------------------------------------------!
  !     Hamiltonian and overlap matrix set up     !
  !-----------------------------------------------!
  CALL timesec(ts0)

  ALLOCATE(h(nmatp,nmatp),o(nmatp,nmatp))

  ! Hamiltonian
  CALL hmlfv(nmatp,ngp,igpig,vgpc,apwalm,h)

  ! overlap
  CALL olpfv(nmatp,ngp,igpig,apwalm,o)

  CALL timesec(ts1)

  timemat = timemat + ts1 - ts0
  
  !---------------------------------------!
  !     solve the eigenvalue equation     !
  !---------------------------------------!
  IF(tefvr) THEN 
    ! system has inversion symmetry: use real symmetric matrix eigen solver
    CALL eveqnfvr(nmatp,ngp,vpc,h,o,evalfv,evecfv)
  ELSE 
    ! no inversion symmetry: use complex Hermitian matrix eigen solver
    CALL eveqnfvz(nmatp,h,o,evalfv,evecfv)
  ENDIF 
  
  DEALLOCATE(h,o)
  
  RETURN 
END SUBROUTINE 

