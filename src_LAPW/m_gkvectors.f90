module m_gkvectors

!------------------------------!
!     G+k-vector variables     !
!------------------------------!
! species for which the muffin-tin radius will be used for calculating gkmax
integer isgkmax
! smallest muffin-tin radius times gkmax
real(8) rgkmax
! maximum |G+k| cut-off for APW functions
real(8) gkmax
! number of G+k-vectors for augmented plane waves
integer, allocatable :: ngk(:,:)
! maximum number of G+k-vectors over all k-points
integer ngkmax
! index from G+k-vectors to G-vectors
integer, allocatable :: igkig(:,:,:)
! G+k-vectors in lattice coordinates
real(8), allocatable :: vgkl(:,:,:,:)
! G+k-vectors in Cartesian coordinates
real(8), allocatable :: vgkc(:,:,:,:)
! length of G+k-vectors
real(8), allocatable :: gkc(:,:,:)
! structure factors for the G+k-vectors
complex(8), allocatable :: sfacgk(:,:,:,:)

end module