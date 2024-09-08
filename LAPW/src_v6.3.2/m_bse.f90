module m_bse

!-------------------------------------------------!
!     Bethe-Salpeter equation (BSE) variables     !
!-------------------------------------------------!
! number of valence and conduction states for transitions
integer nvbse,ncbse
! default number of valence and conduction states
integer nvbse0,ncbse0
! maximum number of extra valence and conduction states
integer, parameter :: maxxbse=20
! number of extra valence and conduction states
integer nvxbse,ncxbse
! extra valence and conduction states
integer istxbse(maxxbse),jstxbse(maxxbse)
! total number of transitions
integer nvcbse
! size of blocks in BSE Hamiltonian matrix
integer nbbse
! size of BSE matrix (= 2*nbbse)
integer nmbse
! index from BSE valence states to second-variational states
integer, allocatable :: istbse(:,:)
! index from BSE conduction states to second-variational states
integer, allocatable :: jstbse(:,:)
! index from BSE valence-conduction pair and k-point to location in BSE matrix
integer, allocatable :: ijkbse(:,:,:)
! BSE Hamiltonian
complex(8), allocatable :: hmlbse(:,:)
! BSE Hamiltonian eigenvalues
real(8), allocatable :: evalbse(:)
! if bsefull is .true. then the full BSE Hamiltonian is calculated, otherwise
! only the Hermitian block
logical bsefull
! if hxbse/hdbse is .true. then the exchange/direct term is included in the BSE
! Hamiltonian
logical hxbse,hdbse

end module