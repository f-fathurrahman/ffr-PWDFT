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
!!   Global variable `Rhoe_R`
!!
SUBROUTINE calc_Rhoe_R( Focc, psi )

  USE m_realspace, ONLY : Npoints, Ns, dVol
  USE m_PWGrid, ONLY : Ngwx, idx_gw2r
  USE m_states, ONLY : Nstates
  USE m_hamiltonian, ONLY : Rhoe_R
  USE m_constants, ONLY : ZZERO
  IMPLICIT NONE
  COMPLEX(8) :: psi(Ngwx,Nstates)
  REAL(8) :: Focc(Nstates)
  COMPLEX(8), ALLOCATABLE :: psiR(:,:)
  INTEGER :: ist, igw, idxG, ip

  ALLOCATE( psiR(Npoints,Nstates) )

  psiR(:,:) = ZZERO

  DO ist = 1, Nstates
    DO igw = 1,Ngwx
      idxG = idx_gw2r(igw)
      psiR(idxG,ist) = psi(igw,ist)
    ENDDO 
    CALL fft_fftw3( psiR(:,ist), Ns(1), Ns(2), Ns(3), .TRUE. )
  ENDDO

  CALL z_ortho_gram_schmidt( psiR, Npoints, Npoints, Nstates )
  psiR(:,:) = psiR(:,:) / sqrt(dVol)

  Rhoe_R(:) = 0.d0
  DO ist = 1,Nstates
    DO ip = 1,Npoints
      Rhoe_R(ip) = Rhoe_R(ip) + Focc(ist) * real( conjg(psiR(ip,ist))*psiR(ip,ist), kind=8)
    ENDDO 
  ENDDO

!  WRITE(*,*)
!  WRITE(*,*) 'Calculating electron density'
!  WRITE(*,'(1x,A,F18.10)') 'Integrated electron density:', sum( Rhoe_R(:) )*dVol

  DEALLOCATE( psiR )

END SUBROUTINE 

