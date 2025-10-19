! to be called from ElkDFTWrapper

!--------------------------------------
subroutine driver_getevecfv(ik, evecfv)
!--------------------------------------
  USE m_states, ONLY: nstfv
  USE m_spin, only: nspnfv
  USE m_kpoints, ONLY: vkl
  USE m_gkvectors, ONLY: vgkl
  USE m_hamiltonian, ONLY: nmatmax
  USE m_misc, ONLY: filext
  !
  implicit none
  integer :: ik
  complex(8) :: evecfv(nmatmax, nstfv, nspnfv)

  ! get the eigenvectors from file
  CALL my_getevecfv(filext, ik, vkl(:,ik), vgkl(:,:,:,ik), evecfv)
  return
end subroutine