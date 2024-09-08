module m_core_states

!------------------------------!
!     core state variables     !
!------------------------------!
! occupancies for core states
REAL(8), ALLOCATABLE :: occcr(:,:)

! eigenvalues for core states
REAL(8), ALLOCATABLE :: evalcr(:,:)

! radial wavefunctions for core states
REAL(8), ALLOCATABLE :: rwfcr(:,:,:,:)

! radial charge density for core states
REAL(8), ALLOCATABLE :: rhocr(:,:,:)

! spincore is .true. if the core is to be treated as spin-polarised
logical spincore

! number of core spin-channels
INTEGER nspncr

END MODULE 
