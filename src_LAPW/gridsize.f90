! !DESCRIPTION:
!   Finds the ${\bf G}$-vector grid which completely contains the vectors with
!   $G<G_{\rm max}$ and is compatible with the FFT routine. The optimal sizes
!   are given by
!   $$ n_i=\frac{G_{\rm max}|{\bf a}_i|}{\pi}+1, $$
!   where ${\bf a}_i$ is the $i$th lattice vector.
SUBROUTINE gridsize(avec,gmaxvr,ngridg,ngtot,intgv)
! !INPUT/OUTPUT PARAMETERS:
!   avec   : lattice vectors (in,real(3,3))
!   gmaxvr : G-vector cut-off (in,real)
!   ngridg : G-vector grid sizes (out,integer(3))
!   ngtot  : total number of G-vectors (out,integer)
!   intgv  : INTEGER grid intervals for each direction (out,integer(2,3))
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(in) :: avec(3,3),gmaxvr
  INTEGER, intent(out) :: ngridg(3),ngtot,intgv(2,3)
  ! local variables
  REAL(8), parameter :: pi=3.1415926535897932385d0
  
  ! find optimal grid size for potential and density
  ngridg(:)=int(gmaxvr*sqrt(avec(1,:)**2+avec(2,:)**2+avec(3,:)**2)/pi)+1
  
  ! find next largest FFT-compatible grid size
  CALL nfftifc(ngridg(1))
  CALL nfftifc(ngridg(2))
  CALL nfftifc(ngridg(3))

  IF((ngridg(1) <= 0).or.(ngridg(2) <= 0).or.(ngridg(3) <= 0)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(gridsize): invalid ngridg : ",3I8)') ngridg
    WRITE(*,*)
    STOP 
  ENDIF 

  ! total number of points in grid
  ngtot = ngridg(1)*ngridg(2)*ngridg(3)

  ! determine INTEGER ranges for grid
  intgv(1,:) = ngridg(:)/2-ngridg(:)+1
  intgv(2,:) = ngridg(:)/2

  RETURN 
END SUBROUTINE 
