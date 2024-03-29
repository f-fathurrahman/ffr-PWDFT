module m_qpoints

!---------------------------!
!     q-point variables     !
!---------------------------!
! q-point grid sizes
integer ngridq(3)
! integer grid intervals for the q-points
integer intq(2,3)
! type of reduction to perform on q-point set (see reducek)
integer reduceq
! number of point group symmetries used for q-point reduction
integer nsymqpt
! point group symmetry matrices used for q-point reduction
integer symqpt(3,3,48)
! total number of reduced q-points
integer nqpt
! total number of non-reduced q-points
integer nqptnr
! map from non-reduced grid to reduced index
integer, allocatable :: iqmap(:,:,:)
! map from non-reduced grid to non-reduced index
integer, allocatable :: iqmapnr(:,:,:)
! locations of q-points on integer grid
integer, allocatable :: ivq(:,:)
! map from (i1,i2,i3) to q-vector index
integer, allocatable :: ivqiq(:,:,:)
! map from q-vector index to complex-complex FFT array
integer, allocatable :: iqfft(:)
! number of complex FFT elements for real-complex transforms
integer nfqrz
! map from q-point index to real-complex FFT index
integer, allocatable :: ifqrz(:)
! map from real-complex FFT index to q-point index
integer, allocatable :: iqrzf(:)
! q-points in lattice coordinates
real(8), allocatable :: vql(:,:)
! q-points in Cartesian coordinates
real(8), allocatable :: vqc(:,:)
! q-point weights
real(8), allocatable :: wqpt(:)
! weight for each non-reduced q-point
real(8) wqptnr
! index of q = 0 point
integer iq0
! regularised Coulomb Green's function in q-space
real(8), allocatable :: gclq(:)

end module