! input:
! - ik: kpoint index
! - evalfv: eigenvalues
! - evecfv: eigenvectors
! Originally from my_eveqn( ik, evalfv, evecfv )

SUBROUTINE debug_hamiltonian( jspn, ik, apwalm, h, o )
  USE m_atoms, ONLY: natmtot
  USE m_gkvectors, ONLY: ngk, ngkmax, igkig, vgkc, gkc, sfacgk
  USE m_hamiltonian, ONLY: nmat, tefvr
  USE m_kpoints, ONLY: vkc
  USE m_muffin_tins, ONLY: lmmaxapw
  USE m_apwlo, ONLY: apwordmax
  USE m_spin, ONLY: nspnfv
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: ik, jspn
  complex(8), intent(out) :: h(nmat(jspn,ik), nmat(jspn,ik))
  complex(8), intent(out) :: o(nmat(jspn,ik), nmat(jspn,ik))
  complex(8), intent(out) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
  !
  ! Local variables
  integer :: ias, i

  write(*,*)
  write(*,*) '<div> ENTER debug_hamiltonian'
  write(*,*)

  ! loop over first-variational spins (nspnfv=2 for spin-spirals only)

  write(*,'(1x,A,I5,A,I5)') 'ik = ', ik, ' jspn = ', jspn
  ! find the matching coefficients
  CALL match( ngk(jspn,ik), vgkc(:,:,jspn,ik), gkc(:,jspn,ik), &
                sfacgk(:,:,jspn,ik), apwalm(:,:,:,:,jspn) )


  ! from hmlfv
  ! zero the upper triangular part of the matrix
  DO i=1,nmat(jspn,ik)
    h(1:i,i) = 0.d0
  ENDDO 
  ! additional zeros to avoid garbage values
  h(:,:) = 0.d0

  ! Debug for 1 atom only
  !ias = 2
  !CALL my_hmlaa(tefvr, ias, ngk(jspn,ik), apwalm(:,:,:,ias,jspn), nmat(jspn,ik), h)
  !
  DO ias=1,natmtot
    CALL hmlaa(tefvr, ias, ngk(jspn,ik), apwalm(:,:,:,ias,jspn), nmat(jspn,ik), h)
  ENDDO 
  CALL my_hmlistl(ngk(jspn,ik), igkig(:,jspn,ik), vgkc(:,:,jspn,ik), nmat(jspn,ik), h)
  !!  
  DO ias = 1,natmtot
    CALL hmlalo(ias, ngk(jspn,ik), apwalm(:,:,:,ias,jspn), nmat(jspn,ik), h)
    CALL hmllolo(ias, ngk(jspn,ik), nmat(jspn,ik), h)
  ENDDO


  ! zero the upper triangular part of the matrix
  DO i = 1,nmat(jspn,ik)
    o(1:i,i) = 0.d0
  ENDDO 
  o(:,:) = 0.d0
  !  
  DO ias=1,natmtot
    CALL olpaa(tefvr, ias, ngk(jspn,ik), apwalm(:,:,:,ias,jspn), nmat(jspn,ik), o)
  ENDDO 
  CALL olpistl(ngk(jspn,ik), igkig(:,jspn,ik), nmat(jspn,ik), o)
  !
  DO ias = 1,natmtot
    CALL olpalo(ias, ngk(jspn,ik), apwalm(:,:,:,ias,jspn), nmat(jspn,ik), o)
    CALL olplolo(ias, ngk(jspn,ik), nmat(jspn,ik), o)
  ENDDO 


  write(*,*)
  write(*,*) '</div> EXIT debug_hamiltonian'
  write(*,*)


  RETURN
END SUBROUTINE 
