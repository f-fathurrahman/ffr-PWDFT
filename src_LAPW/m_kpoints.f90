module m_kpoints

!---------------------------!
!     k-point variables     !
!---------------------------!
! autokpt is .true. if the k-point set is determined automatically
logical autokpt
! radius of sphere used to determine k-point density when autokpt is .true.
real(8) radkpt
! k-point grid sizes
integer ngridk(3)
! k-point offset
real(8) vkloff(3)
! corners of box in lattice coordinates containing the k-points, the zeroth
! vector is the origin
real(8) kptboxl(3,0:3)
! type of reduction to perform on k-point set
!  0 : no reduction
!  1 : reduce with full crystal symmetry group
!  2 : reduce with symmorphic symmetries only
integer reducek
! number of point group symmetries used for k-point reduction
integer nsymkpt
! point group symmetry matrices used for k-point reduction
integer symkpt(3,3,48)
! total number of reduced k-points
integer nkpt
! total number of non-reduced k-points
integer nkptnr
! locations of k-points on integer grid
integer, allocatable :: ivk(:,:)
! map from (i1,i2,i3) to reduced k-point index
integer, allocatable :: ivkik(:,:,:)
! map from (i1,i2,i3) to non-reduced k-point index
integer, allocatable :: ivkiknr(:,:,:)
! k-points in lattice coordinates
real(8), allocatable :: vkl(:,:)
! k-points in Cartesian coordinates
real(8), allocatable :: vkc(:,:)
! reduced k-point weights
real(8), allocatable :: wkpt(:)
! weight of each non-reduced k-point
real(8) wkptnr
! k-point at which to determine effective mass tensor
real(8) vklem(3)
! displacement size for computing the effective mass tensor
real(8) deltaem
! number of displacements in each direction
integer ndspem

end module 
