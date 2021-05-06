module m_muffin_tins

use m_atoms, only: maxspecies

!---------------------------------------------------------------!
!     muffin-tin radial mesh and angular momentum variables     !
!---------------------------------------------------------------!
! scale factor for number of muffin-tin points
real(8) nrmtscf
! number of muffin-tin radial points for each species
integer nrmt(maxspecies)
! maximum nrmt over all the species
integer nrmtmax
! optional default muffin-tin radius for all atoms
real(8) rmtall
! minimum allowed distance between muffin-tin surfaces
real(8) rmtdelta
! muffin-tin radii
real(8) rmt(maxspecies)
! total muffin-tin volume
real(8) omegamt
! radial step length for coarse mesh
integer lradstp
! number of coarse radial mesh points
integer nrcmt(maxspecies)
! maximum nrcmt over all the species
integer nrcmtmax
! coarse muffin-tin radial mesh
real(8), allocatable :: rcmt(:,:)
! r^l on fine radial mesh
real(8), allocatable :: rlmt(:,:,:)
! r^l on coarse radial mesh
real(8), allocatable :: rlcmt(:,:,:)
! weights for spline integration on fine radial mesh
real(8), allocatable :: wrmt(:,:)
! weights for spline partial integration on fine radial mesh
real(8), allocatable :: wprmt(:,:,:)
! weights for spline integration on coarse radial mesh
real(8), allocatable :: wrcmt(:,:)
! weights for spline partial integration on coarse radial mesh
real(8), allocatable :: wprcmt(:,:,:)
! maximum allowable angular momentum for augmented plane waves
integer, parameter :: maxlapw=50
! maximum angular momentum for augmented plane waves
integer lmaxapw
! (lmaxapw+1)^2
integer lmmaxapw
! maximum angular momentum on the outer part of the muffin-tin
integer lmaxo
! (lmaxo+1)^2
integer lmmaxo
! maximum angular momentum on the inner part of the muffin-tin
integer lmaxi
! (lmaxi+1)^2
integer lmmaxi
! fraction of muffin-tin radius which constitutes the inner part
real(8) fracinr
! number of fine/coarse radial points on the inner part of the muffin-tin
integer nrmti(maxspecies),nrcmti(maxspecies)
! index to (l,m) pairs
integer, allocatable :: idxlm(:,:)
! inverse index to (l,m) pairs
integer, allocatable :: idxil(:),idxim(:)
! number of fine/coarse points in packed muffin-tins
integer npmti(maxspecies),npmt(maxspecies)
integer npcmti(maxspecies),npcmt(maxspecies)
! maximum number of points over all packed muffin-tins
integer npmtmax,npcmtmax


end module