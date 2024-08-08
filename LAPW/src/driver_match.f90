!----------------------------------------
subroutine driver_match(ik, apwalm)
!----------------------------------------
  USE m_atoms, ONLY: natmtot
  USE m_gkvectors, ONLY: ngk, ngkmax, vgkc, gkc, sfacgk
  USE m_muffin_tins, ONLY: lmmaxapw
  USE m_apwlo, ONLY: apwordmax
  USE m_spin, only: nspnfv
  implicit none
  ! argument
  integer :: ik
  COMPLEX(8) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
  ! Local variables
  integer :: jspn

  write(*,*)
  write(*,*) '>>>>>    ENTER driver_match'
  write(*,*)

  ! find the matching coefficients
  do jspn = 1,nspnfv
    CALL match( ngk(jspn,ik), vgkc(:,:,jspn,ik), gkc(:,jspn,ik), &
                sfacgk(:,:,jspn,ik), apwalm(:,:,:,:,jspn) )
    write(*,*) 'jspn = ', jspn, ' is done'
  enddo

  write(*,*)
  write(*,*) '>>>>>    EXIT driver_match'
  write(*,*)

  return
end subroutine
