module m_hamiltonian

!-------------------------------------------!
!     overlap and Hamiltonian variables     !
!-------------------------------------------!
! overlap and Hamiltonian matrices sizes at each k-point
integer, allocatable :: nmat(:,:)
! maximum nmat over all k-points
integer nmatmax
! index to the position of the local-orbitals in the H and O matrices
integer, allocatable :: idxlo(:,:,:)
! APW-local-orbital overlap integrals
real(8), allocatable :: oalo(:,:,:)
! local-orbital-local-orbital overlap integrals
real(8), allocatable :: ololo(:,:,:)
! APW-APW Hamiltonian integrals
real(8), allocatable :: haa(:,:,:,:,:,:)
! local-orbital-APW Hamiltonian integrals
real(8), allocatable :: hloa(:,:,:,:,:)
! local-orbital-local-orbital Hamiltonian integrals
real(8), allocatable :: hlolo(:,:,:,:)
! complex Gaunt coefficient array
complex(8), allocatable :: gntyry(:,:,:)
! tefvr is .true. if the first-variational eigenvalue equation is to be solved
! as a real symmetric problem
logical tefvr
! tefvit is .true. if the first-variational eigenvalue equation is to be solved
! iteratively
logical tefvit
! minimum and maximum allowed number of eigenvalue equation iterations
integer minitefv,maxitefv
! eigenvalue mixing parameter for iterative solver
real(8) befvit
! iterative solver convergence tolerance
real(8) epsefvit
! type of eigenvalue solver to be used
integer evtype

end module
