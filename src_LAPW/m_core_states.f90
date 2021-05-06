module m_core_states

!------------------------------!
!     core state variables     !
!------------------------------!
! occupancies for core states
real(8), allocatable :: occcr(:,:)
! eigenvalues for core states
real(8), allocatable :: evalcr(:,:)
! radial wavefunctions for core states
real(8), allocatable :: rwfcr(:,:,:,:)
! radial charge density for core states
real(8), allocatable :: rhocr(:,:,:)
! spincore is .true. if the core is to be treated as spin-polarised
logical spincore
! number of core spin-channels
integer nspncr

end module