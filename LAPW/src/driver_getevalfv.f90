! to be called from ElkDFTWrapper

!--------------------------------------
subroutine driver_getevalfv(ik, evalfv)
!--------------------------------------
  USE m_states, ONLY: nstfv
  USE m_spin, only: nspnfv
  USE m_kpoints, ONLY: vkl
  USE m_hamiltonian, ONLY: nmatmax
  USE m_misc, ONLY: filext
  !
  implicit none
  integer :: ik
  real(8) :: evalfv(nstfv,nspnfv)

  ! get the eigenvectors from file
  CALL getevalfv(filext, ik, vkl(:,ik), evalfv)
  return
end subroutine
