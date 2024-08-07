subroutine driver_match(jspn, ik, apwalm)
  USE m_atoms, ONLY: natmtot
  USE m_gkvectors, ONLY: ngk, ngkmax, vgkc, gkc, sfacgk
  USE m_muffin_tins, ONLY: lmmaxapw
  USE m_apwlo, ONLY: apwordmax
  USE m_spin, only: nspnfv
  implicit none
  ! argument
  integer :: jspn, ik
  COMPLEX(8) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
  
  ! find the matching coefficients
  CALL match( ngk(jspn,ik), vgkc(:,:,jspn,ik), gkc(:,jspn,ik), &
              sfacgk(:,:,jspn,ik), apwalm(:,:,:,:,jspn) )
  return
end subroutine
