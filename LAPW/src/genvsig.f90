SUBROUTINE genvsig()
  USE m_gvectors, ONLY: ngtot, ngridg, cfunir, igfft, ngvec
  USE m_density_pot_xc, ONLY: vsig, vsir
  !
  ! Generates the Fourier transform of the Kohn-Sham effective potential in the
  ! interstitial region. The potential is first multiplied by the characteristic
  ! function which zeros it in the muffin-tins. See routine {\tt gencfun}.
  !
  IMPLICIT NONE 
  ! allocatable arrays
  COMPLEX(8), ALLOCATABLE :: zfft(:)
  ALLOCATE( zfft(ngtot) )
  
  ! multiply potential by characteristic function in real-space
  zfft(:) = vsir(:)*cfunir(:)
  
  ! Fourier transform to G-space
  CALL zfftifc(3, ngridg, -1, zfft)
  
  ! store in global array
  vsig(1:ngvec) = zfft(igfft(1:ngvec))

  DEALLOCATE( zfft )
  RETURN 
END SUBROUTINE 

