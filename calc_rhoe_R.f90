!!
!! PURPOSE:
!!
!!   This subroutine calculates electronic density, given `psi`
!!   (which need not to be Kohn-Sham states) and occupation number `Focc`
!! 
!! AUTHOR:
!!
!!   Fadjar Fathurrahman
!!
!! MODIFY:
!!
!!   Global variable `rhoe`
!!
SUBROUTINE calc_rhoe_R( psiR, Focc )

  USE m_realspace, ONLY : Npoints
  USE m_states, ONLY : Nstates
  USE m_hamiltonian, ONLY : Rhoe
  IMPLICIT NONE
  REAL(8) :: psiR(Npoints,Nstates)
  REAL(8) :: Focc(Nstates)
  INTEGER :: ist

  Rhoe(:) = 0.d0
  DO ist = 1, Nstates
    Rhoe(:) = Rhoe(:) + Focc(ist) * psiR(:,ist) * psiR(:,ist)
  ENDDO

  WRITE(*,*)
  WRITE(*,*) 'Calculating electron density'
  WRITE(*,'(1x,A,F18.10)') 'Integrated electron density:', sum( Rhoe(:) )*dVol

END SUBROUTINE 

