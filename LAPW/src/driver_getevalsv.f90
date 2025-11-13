! to be called from ElkDFTWrapper

!--------------------------------------
subroutine driver_getevalsv(ik, evalsv)
!--------------------------------------
  USE m_states, ONLY: nstsv
  USE m_spin, only: nspnfv
  USE m_kpoints, ONLY: vkl
  USE m_misc, ONLY: filext
  !
  implicit none
  integer :: ik
  real(8) :: evalsv(nstsv,nspnfv)

  ! get the eigenvectors from file
  CALL getevalsv(filext, ik, vkl(:,ik), evalsv)
  return
end subroutine
