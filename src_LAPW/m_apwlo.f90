module m_apwlo

use m_atoms, only: maxspecies
use m_muffin_tins, only: maxlapw

!-----------------------------------------!
!     APW and local-orbital variables     !
!-----------------------------------------!
! energy step used for numerical calculation of energy derivatives
real(8) :: deapwlo
! maximum allowable APW order
integer, parameter :: maxapword=4
! APW order
integer ::apword(0:maxlapw,maxspecies)
! maximum of apword over all angular momenta and species
integer apwordmax
! total number of APW coefficients (l, m and order) for each species
integer lmoapw(maxspecies)
! polynomial order used for APW radial derivatives
integer npapw
! APW initial linearisation energies
real(8) apwe0(maxapword,0:maxlapw,maxspecies)
! APW linearisation energies
real(8), allocatable :: apwe(:,:,:)
! APW derivative order
integer apwdm(maxapword,0:maxlapw,maxspecies)
! apwve is .true. if the linearisation energies are allowed to vary
logical apwve(maxapword,0:maxlapw,maxspecies)
! APW radial functions
real(8), allocatable :: apwfr(:,:,:,:,:)
! derivate of radial functions at the muffin-tin surface
real(8), allocatable :: apwdfr(:,:,:)
! maximum number of local-orbitals
integer, parameter :: maxlorb=200
! maximum allowable local-orbital order
integer, parameter :: maxlorbord=5
! number of local-orbitals
integer nlorb(maxspecies)
! maximum nlorb over all species
integer nlomax
! total number of local-orbitals
integer nlotot
! local-orbital order
integer lorbord(maxlorb,maxspecies)
! maximum lorbord over all species
integer lorbordmax
! polynomial order used for local-orbital radial derivatives
integer nplorb
! local-orbital angular momentum
integer lorbl(maxlorb,maxspecies)
! maximum lorbl over all species
integer lolmax
! (lolmax+1)^2
integer lolmmax
! local-orbital initial energies
real(8) lorbe0(maxlorbord,maxlorb,maxspecies)
! local-orbital energies
real(8), allocatable :: lorbe(:,:,:)
! local-orbital derivative order
integer lorbdm(maxlorbord,maxlorb,maxspecies)
! lorbve is .true. if the linearisation energies are allowed to vary
logical lorbve(maxlorbord,maxlorb,maxspecies)
! local-orbital radial functions
real(8), allocatable :: lofr(:,:,:,:)
! band energy search tolerance
real(8) epsband
! maximum allowed change in energy during band energy search; enforced only if
! default energy is less than zero
real(8) demaxbnd
! minimum default linearisation energy over all APWs and local-orbitals
real(8) e0min
! if autolinengy is .true. then the fixed linearisation energies are set to the
! Fermi energy minus dlefe
logical autolinengy
! difference between linearisation and Fermi energies when autolinengy is .true.
real(8) dlefe
! lorbcnd is .true. if conduction state local-orbitals should be added
logical lorbcnd
! conduction state local-orbital order
integer lorbordc
! excess order of the APW and local-orbital functions
integer nxoapwlo
! excess local orbitals
integer nxlo

end module
