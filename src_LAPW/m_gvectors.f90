module m_gvectors

!----------------------------!
!     G-vector variables     !
!----------------------------!
! G-vector cut-off for interstitial potential and density
real(8) gmaxvr
! G-vector grid sizes
integer ngridg(3)
! G-vector grid sizes for coarse grid (G < 2*gkmax)
integer ngdc(3)
! total number of G-vectors
integer ngtot
! total number of G-vectors for coarse grid (G < 2*gkmax)
integer ngtc
! integer grid intervals for each direction
integer intgv(2,3)
! number of G-vectors with G < gmaxvr
integer ngvec
! number of G-vectors for coarse grid (G < 2*gkmax)
integer ngvc
! G-vector integer coordinates (i1,i2,i3)
integer, allocatable :: ivg(:,:)
! map from (i1,i2,i3) to G-vector index
integer, allocatable :: ivgig(:,:,:)
! map from G-vector index to FFT array
integer, allocatable :: igfft(:)
! map from G-vector index to FFT array for coarse grid (G < 2*gkmax)
integer, allocatable :: igfc(:)
! G-vectors in Cartesian coordinates
real(8), allocatable :: vgc(:,:)
! length of G-vectors
real(8), allocatable :: gc(:)
! Coulomb Green's function in G-space = 4 pi / G^2
real(8), allocatable :: gclg(:)
! spherical Bessel functions j_l(|G|R_mt)
real(8), allocatable :: jlgrmt(:,:,:)
! spherical harmonics of the G-vectors
complex(8), allocatable :: ylmg(:,:)
! structure factors for the G-vectors
complex(8), allocatable :: sfacg(:,:)
! smooth step function form factors for all species and G-vectors
real(8), allocatable :: ffacg(:,:)
! characteristic function in G-space: 0 inside the muffin-tins and 1 outside
complex(8), allocatable :: cfunig(:)
! characteristic function in real-space: 0 inside the muffin-tins and 1 outside
real(8), allocatable :: cfunir(:)
! characteristic function in real-space for coarse grid (G < 2*gkmax)
real(8), allocatable :: cfrc(:)

end module
