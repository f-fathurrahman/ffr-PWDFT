PROGRAM main
  IMPLICIT NONE 

  ! read input and initialize some variables
  CALL read_input()
  CALL init0()
  CALL init1()
  
  CALL info_crystal()
  CALL info_symmetry()
  CALL writesym()
  CALL info_gvectors()
  CALL info_muffin_tins()

  !CALL rhoinit()
  CALL gndstate()

  CALL test_rfplot()

END PROGRAM

!-----------------------
subroutine test_rfplot()
!-----------------------
  USE m_density_pot_xc, ONLY: rhoir, rhomt
  USE m_atoms, ONLY: natmtot, idxis
  USE m_gvectors, ONLY: ngtot, ngridg
  USE m_muffin_tins, ONLY: lmmaxo, nrmt, nrmti, nrmtmax
  implicit none
  REAL(8), ALLOCATABLE :: rfmt1(:,:,:), vpl(:,:), fp(:)
  COMPLEX(8), ALLOCATABLE :: zfft(:)
  integer :: np
  integer :: ias, is, ip

  ! unpack the muffin-tin function
  ALLOCATE(rfmt1(lmmaxo,nrmtmax,natmtot))
  DO ias=1,natmtot
    is=idxis(ias)
    CALL rfmtpack(.false., nrmt(is), nrmti(is), rhomt(:,ias), rfmt1(:,:,ias))
  ENDDO 

  ! Fourier transform rfir to G-space
  ALLOCATE(zfft(ngtot))
  zfft(:) = rhoir(:)
  CALL zfftifc(3,ngridg,-1,zfft)

  np = 2
  allocate( vpl(3,np) )
  allocate( fp(np) )
  vpl(1:3,1) = (/ 0.0d0, 0.0d0, 0.01d0 /)
  vpl(1:3,2) = (/ 0.25d0, 0.25d0, 0.20d0 /)

  ! begin loop over all points
  DO ip=1,np
    CALL my_rfip(ip, np, vpl, rfmt1, zfft, fp)
    write(*,*) 'fp = ', fp(ip)
  ENDDO 

  DEALLOCATE(rfmt1,zfft)
  deallocate(vpl, fp)

  RETURN 
end subroutine

include 'my_rfip.f90'