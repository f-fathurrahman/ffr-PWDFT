module m_states

!--------------------------------------------!
!     eigenvalue and occupancy variables     !
!--------------------------------------------!
! number of empty states per atom and spin
real(8) nempty0
! number of empty states
integer nempty
! number of first-variational states
integer nstfv
! number of second-variational states
integer nstsv
! smearing type
integer stype
! smearing function description
character(256) sdescr
! smearing width
real(8) swidth
! autoswidth is .true. if the smearing width is to be determined automatically
logical autoswidth
! effective mass used in smearing width formula
real(8) mstar
! maximum allowed occupancy (1 or 2)
real(8) occmax
! convergence tolerance for occupancies
real(8) epsocc
! second-variational occupation numbers
real(8), allocatable :: occsv(:,:)
! Fermi energy for second-variational states
real(8) efermi
! scissor correction applied when computing response functions
real(8) scissor
! density of states at the Fermi energy
real(8) fermidos
! estimated indirect and direct band gaps
real(8) bandgap(2)
! k-points of indirect and direct gaps
integer ikgap(3)
! error tolerance for the first-variational eigenvalues
real(8) evaltol
! second-variational eigenvalues
real(8), allocatable :: evalsv(:,:)
! tevecsv is .true. if second-variational eigenvectors are calculated
logical tevecsv
! maximum number of k-point and states indices in user-defined list
integer, parameter :: maxkst=20
! number of k-point and states indices in user-defined list
integer nkstlist
! user-defined list of k-point and state indices
integer kstlist(2,maxkst)

end module