module m_oep_hf

!----------------------------------------!
!     OEP and Hartree-Fock variables     !
!----------------------------------------!
! maximum number of core states over all species
integer ncrmax
! maximum number of OEP iterations
integer maxitoep
! OEP step size
real(8) tauoep
! magnitude of the OEP residual
real(8) resoep
! exchange potential and magnetic field
real(8), allocatable :: vxmt(:,:),vxir(:)
real(8), allocatable :: bxmt(:,:,:),bxir(:,:)
! hybrid is .true. if a hybrid functional is to be used
logical hybrid
! hybrid functional mixing coefficient
real(8) hybridc

end module