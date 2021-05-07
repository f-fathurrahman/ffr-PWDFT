SUBROUTINE symrvf(tspin,tnc,nr,nri,np,ld,rvfmt,rvfir)
  ! !USES:
  USE modmain
  ! !INPUT/OUTPUT PARAMETERS:
  !   tspin : .true. if spin rotations should be used (in,logical)
  !   tnc   : .true. if the vector field is non-collinear, otherwise it is
  !           collinear along the z-axis (in,logical)
  !   nr    : number of radial points for each species (in,integer(nspecies))
  !   nri   : number of radial points on the inner part (in,integer(nspecies))
  !   np    : total number of points in each muffin-tin (in,integer(nspecies))
  !   ld    : leading dimension (in,integer)
  !   rvfmt : real muffin-tin vector field (in,real(ld,natmtot,*))
  !   rvfir : real interstitial vector field (in,real(ngtot,*))
  ! !DESCRIPTION:
  !   Symmetrises a vector field defined over the entire unit cell using the full
  !   set of crystal symmetries. If a particular symmetry involves rotating atom
  !   1 into atom 2, THEN  the spatial and spin rotations of that symmetry are
  !   applied to the vector field in atom 2 (expressed in spherical harmonic
  !   coefficients), which is THEN  added to the field in atom 1. This is repeated
  !   for all symmetry operations. The fully symmetrised field in atom 1 is THEN 
  !   rotated and copied to atom 2. Symmetrisation of the interstitial part of the
  !   field is performed by {\tt symrvfir}. See also {\tt symrfmt} and
  !   {\tt findsym}.

  IMPLICIT NONE 
  ! arguments
  logical, intent(in) :: tspin,tnc
  INTEGER, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
  INTEGER, intent(in) :: ld
  REAL(8), intent(inout) :: rvfmt(ld,natmtot,*),rvfir(ngtot,*)
  ! local variables
  CALL symrvfmt(tspin,tnc,nr,nri,np,ld,rvfmt)
  CALL symrvfir(tspin,tnc,rvfir)
  RETURN 
END SUBROUTINE 
