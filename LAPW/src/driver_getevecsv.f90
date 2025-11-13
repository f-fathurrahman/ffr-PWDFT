! to be called from ElkDFTWrapper

!--------------------------------------
subroutine driver_getevecsv(ik, evecsv)
!--------------------------------------
  USE m_states, ONLY: nstsv
  USE m_kpoints, ONLY: vkl
  USE m_misc, ONLY: filext
  !
  implicit none
  integer :: ik
  complex(8) :: evecsv(nstsv,nstsv)

  ! get the eigenvectors from file
  CALL getevecsv(filext, ik, vkl(:,ik), evecsv)
  return
end subroutine