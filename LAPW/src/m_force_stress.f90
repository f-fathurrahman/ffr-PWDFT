module m_force_stress

!------------------------------------!
!     force and stress variables     !
!------------------------------------!
! tforce is .true. if force should be calculated
logical tforce
! Hellmann-Feynman force on each atom
real(8), allocatable :: forcehf(:,:)
! incomplete basis set (IBS) force on each atom
real(8), allocatable :: forceibs(:,:)
! total force on each atom
real(8), allocatable :: forcetot(:,:)
! previous total force on each atom
real(8), allocatable :: forcetotp(:,:)
! maximum force magnitude over all atoms
real(8) forcemax
! tfav0 is .true. if the average force should be zero in order to prevent
! translation of the atomic basis
logical tfav0
! average force
real(8) forceav(3)
! atomic position optimisation type
!  0 : no optimisation
!  1 : unconstrained optimisation
integer atpopt
! maximum number of atomic position optimisation steps
integer maxatpstp
! default step size parameter for atomic position optimisation
real(8) tau0atp
! step size parameters for each atom
real(8), allocatable :: tauatp(:)
! number of strain tensors
integer nstrain
! current strain tensor
integer istrain
data istrain / 0 /
! strain tensors
real(8) strain(3,3,9)
! infinitesimal displacement parameter multiplied by the strain tensor for
! computing the stress tensor
real(8) deltast
! symmetry reduced stress tensor components
real(8) stress(9)
! previous stress tensor
real(8) stressp(9)
! stress tensor component magnitude maximum
real(8) stressmax
! lattice vector optimisation type
!  0 : no optimisation
!  1 : unconstrained optimisation
!  2 : iso-volumetric optimisation
integer latvopt
! maximum number of lattice vector optimisation steps
integer maxlatvstp
! default step size parameter for lattice vector optimisation
real(8) tau0latv
! step size for each stress tensor component acting on the lattice vectors
real(8) taulatv(9)

end module
