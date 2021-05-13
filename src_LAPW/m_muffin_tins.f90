module m_muffin_tins

use m_atoms, only: maxspecies

!---------------------------------------------------------------!
!     muffin-tin radial mesh and angular momentum variables     !
!---------------------------------------------------------------!
! scale factor for number of muffin-tin points
REAL(8) :: nrmtscf
! number of muffin-tin radial points for each species
INTEGER :: nrmt(maxspecies)
! maximum nrmt over all the species
INTEGER :: nrmtmax
! optional default muffin-tin radius for all atoms
REAL(8) :: rmtall
! minimum allowed distance between muffin-tin surfaces
REAL(8) :: rmtdelta
! muffin-tin radii
REAL(8) :: rmt(maxspecies)
! total muffin-tin volume
REAL(8) :: omegamt
! radial step length for coarse mesh
INTEGER :: lradstp
! number of coarse radial mesh points
INTEGER :: nrcmt(maxspecies)
! maximum nrcmt over all the species
INTEGER :: nrcmtmax
! coarse muffin-tin radial mesh
REAL(8), ALLOCATABLE :: rcmt(:,:)
! r^l on fine radial mesh
REAL(8), ALLOCATABLE :: rlmt(:,:,:)
! r^l on coarse radial mesh
REAL(8), ALLOCATABLE :: rlcmt(:,:,:)
! weights for spline integration on fine radial mesh
REAL(8), ALLOCATABLE :: wrmt(:,:)
! weights for spline partial integration on fine radial mesh
REAL(8), ALLOCATABLE :: wprmt(:,:,:)
! weights for spline integration on coarse radial mesh
REAL(8), ALLOCATABLE :: wrcmt(:,:)
! weights for spline partial integration on coarse radial mesh
REAL(8), ALLOCATABLE :: wprcmt(:,:,:)
! maximum allowable angular momentum for augmented plane waves
INTEGER, parameter :: maxlapw=50
! maximum angular momentum for augmented plane waves
INTEGER lmaxapw
! (lmaxapw+1)^2
INTEGER lmmaxapw
! maximum angular momentum on the outer part of the muffin-tin
INTEGER lmaxo
! (lmaxo+1)^2
INTEGER lmmaxo
! maximum angular momentum on the inner part of the muffin-tin
INTEGER lmaxi
! (lmaxi+1)^2
INTEGER lmmaxi
! fraction of muffin-tin radius which constitutes the inner part
REAL(8) fracinr
! number of fine/coarse radial points on the inner part of the muffin-tin
INTEGER nrmti(maxspecies),nrcmti(maxspecies)
! index to (l,m) pairs
INTEGER, ALLOCATABLE :: idxlm(:,:)
! inverse index to (l,m) pairs
INTEGER, ALLOCATABLE :: idxil(:),idxim(:)
! number of fine/coarse points in packed muffin-tins
INTEGER npmti(maxspecies),npmt(maxspecies)
INTEGER npcmti(maxspecies),npcmt(maxspecies)
! maximum number of points over all packed muffin-tins
INTEGER npmtmax,npcmtmax


end module