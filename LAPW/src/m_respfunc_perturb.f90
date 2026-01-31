module m_respfunc_perturb

implicit none

!-------------------------------------------------------------!
!     response function and perturbation theory variables     !
!-------------------------------------------------------------!
! |G| cut-off for response functions
real(8) gmaxrf
! energy cut-off for response functions
real(8) emaxrf
! number of G-vectors for response functions
integer ngrf
! number of response function frequencies
integer nwrf
! complex response function frequencies
complex(8), allocatable :: wrf(:)
! maximum number of spherical Bessel functions on the coarse radial mesh over
! all species
integer njcmax

end module
