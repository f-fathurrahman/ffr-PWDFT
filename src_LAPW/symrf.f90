SUBROUTINE symrf(nr,nri,np,ld,rfmt,rfir)
  ! !USES:
  USE m_atoms, ONLY: nspecies, natmtot
  USE m_gvectors, ONLY: ngtot
  ! !INPUT/OUTPUT PARAMETERS:
  !   nr   : number of radial points for each species (in,integer(nspecies))
  !   nri  : number of radial points on the inner part (in,integer(nspecies))
  !   np   : total number of points in each muffin-tin (in,integer(nspecies))
  !   ld   : leading dimension (in,integer)
  !   rfmt : real muffin-tin function (inout,real(ld,natmtot))
  !   rfir : real intersitial function (inout,real(ngtot))
  ! !DESCRIPTION:
  !   Symmetrises a real scalar function defined over the entire unit cell using
  !   the full set of crystal symmetries. In the muffin-tin of a particular atom
  !   the spherical harmonic coefficients of every equivlent atom are rotated and
  !   averaged. The interstitial part of the function is first Fourier transformed
  !   to $G$-space, and THEN  averaged over each symmetry by rotating the Fourier
  !   coefficients and multiplying them by a phase factor corresponding to the
  !   symmetry translation.
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
  INTEGER, intent(in) :: ld
  REAL(8), intent(inout) :: rfmt(ld,natmtot),rfir(ngtot)
  ! local variables
  CALL symrfmt(nr,nri,np,ld,rfmt)
  CALL symrfir(rfir)
  RETURN 
END SUBROUTINE 
