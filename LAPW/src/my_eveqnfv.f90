!-------------------------------------------------------------------
SUBROUTINE my_eveqnfv(nmatp,ngp,igpig,vpc,vgpc,apwalm,evalfv,evecfv)
!-------------------------------------------------------------------
  USE m_atoms, ONLY: natmtot
  USE m_gkvectors, ONLY: ngkmax
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
  ! allocatable arrays
  COMPLEX(8), allocatable :: h(:,:),o(:,:)


  write(*,*)
  write(*,*) '<div> ENTER my_eveqnfv'
  write(*,*)

  write(*,*) 'nmatp = ', nmatp ! may include local orbital
  write(*,*) 'ngp = ', ngp ! should be the same as nmatp in case no local orbitals are used
  ! igpig, vpc, vgpc: igk index to ig, k-vector, G+k-vector

  !-----------------------------------------------!
  !     Hamiltonian and overlap matrix set up     !
  !-----------------------------------------------!
  ALLOCATE(h(nmatp,nmatp),o(nmatp,nmatp))

  ! Hamiltonian
  CALL hmlfv(nmatp,ngp,igpig,vgpc,apwalm,h)

  ! overlap
  CALL olpfv(nmatp,ngp,igpig,apwalm,o)

  
  !---------------------------------------!
  !     solve the eigenvalue equation     !
  !---------------------------------------!
  IF(tefvr) THEN 
    write(*,*) '*** Using real symmetric matrix eigensolver'
    ! system has inversion symmetry: use real symmetric matrix eigen solver
    CALL eveqnfvr(nmatp,ngp,vpc,h,o,evalfv,evecfv)
  ELSE 
    write(*,*) '*** Using complex Hermitian matrix eigensolver'
    ! no inversion symmetry: use complex Hermitian matrix eigen solver
    CALL eveqnfvz(nmatp,h,o,evalfv,evecfv)
  ENDIF 
  
  DEALLOCATE(h,o)

  write(*,*)
  write(*,*) '</div> EXIT my_eveqnfv'
  write(*,*)

  
  RETURN 
END SUBROUTINE 

