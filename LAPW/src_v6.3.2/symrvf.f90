
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrvf
! !INTERFACE:
subroutine symrvf(tspin,tnc,nr,nri,np,ld,rvfmt,rvfir)
! !USES:
use modmain
use modomp
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
!   1 into atom 2, then the spatial and spin rotations of that symmetry are
!   applied to the vector field in atom 2 (expressed in spherical harmonic
!   coefficients), which is then added to the field in atom 1. This is repeated
!   for all symmetry operations. The fully symmetrised field in atom 1 is then
!   rotated and copied to atom 2. Symmetrisation of the interstitial part of the
!   field is performed by {\tt symrvfir}. See also {\tt symrfmt} and
!   {\tt findsym}.
!
! !REVISION HISTORY:
!   Created May 2007 (JKD)
!   Fixed problem with improper rotations, February 2008 (L. Nordstrom,
!    F. Bultmark and F. Cricchio)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tspin,tnc
integer, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
integer, intent(in) :: ld
real(8), intent(inout) :: rvfmt(ld,natmtot,*),rvfir(ngtot,*)
! local variables
integer nthd
call holdthd(2,nthd)
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP SECTION
call symrvfmt(tspin,tnc,nr,nri,np,ld,rvfmt)
!$OMP SECTION
call symrvfir(tspin,tnc,rvfir)
!$OMP END PARALLEL SECTIONS
call freethd(nthd)
return
end subroutine
!EOC

