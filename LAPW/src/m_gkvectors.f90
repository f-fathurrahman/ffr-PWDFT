MODULE m_gkvectors

!------------------------------!
!     G+k-vector variables     !
!------------------------------!
! species for which the muffin-tin radius will be used for calculating gkmax
INTEGER isgkmax
! smallest muffin-tin radius times gkmax
REAL(8) rgkmax
! maximum |G+k| cut-off for APW functions
REAL(8) gkmax
! number of G+k-vectors for augmented plane waves
INTEGER, ALLOCATABLE :: ngk(:,:)
! maximum number of G+k-vectors over all k-points
INTEGER ngkmax
! index from G+k-vectors to G-vectors
INTEGER, ALLOCATABLE :: igkig(:,:,:)
! G+k-vectors in lattice coordinates
REAL(8), ALLOCATABLE :: vgkl(:,:,:,:)
! G+k-vectors in Cartesian coordinates
REAL(8), ALLOCATABLE :: vgkc(:,:,:,:)
! length of G+k-vectors
REAL(8), ALLOCATABLE :: gkc(:,:,:)
! structure factors for the G+k-vectors
COMPLEX(8), ALLOCATABLE :: sfacgk(:,:,:,:)

END MODULE 