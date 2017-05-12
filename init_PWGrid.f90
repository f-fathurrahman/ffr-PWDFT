SUBROUTINE init_PWGrid( ecutwfc_in, LatVecs_in )
  USE m_constants, ONLY : PI
  USE m_realspace, ONLY : Ns, Npoints, dVol
  USE m_cell, ONLY : LatVecs, RecVecs, CellVolume, LatVecsLen
  USE m_PWGrid, ONLY : ecutwfc, ecutrho
  USE fft_support, ONLY : good_fft_order 
  IMPLICIT NONE 
  REAL(8) :: ecutwfc_in
  REAL(8) :: LatVecs_in(3,3)
  ! function
  REAL(8) :: det_m3x3

  ecutwfc = ecutwfc_in
  ecutrho = 4.d0*ecutwfc

  LatVecs(:,:) = LatVecs_in(:,:)
  CALL inv_m3x3( transpose(LatVecs), RecVecs )
  RecVecs(:,:) = 2.d0*PI*RecVecs(:,:)
  CellVolume = det_m3x3( LatVecs )

  LatVecsLen(1) = sqrt( LatVecs(1,1)**2 + LatVecs(1,2)**2 + LatVecs(1,3)**2 )
  LatVecsLen(2) = sqrt( LatVecs(2,1)**2 + LatVecs(2,2)**2 + LatVecs(2,3)**2 )
  LatVecsLen(3) = sqrt( LatVecs(3,1)**2 + LatVecs(3,2)**2 + LatVecs(3,3)**2 )

  Ns(1) = good_fft_order( 2*nint( sqrt(0.5d0*ecutrho)*LatVecsLen(1)/PI ) + 1 )
  Ns(2) = good_fft_order( 2*nint( sqrt(0.5d0*ecutrho)*LatVecsLen(2)/PI ) + 1 )
  Ns(3) = good_fft_order( 2*nint( sqrt(0.5d0*ecutrho)*LatVecsLen(3)/PI ) + 1 )

  Npoints = Ns(1) * Ns(2) * Ns(3)

  CALL init_gvec()
  CALL init_gvecw()

  dVol = CellVolume/Npoints

END SUBROUTINE 

