module m_apwlo

use m_atoms, only: maxspecies
use m_muffin_tins, only: maxlapw

!-----------------------------------------!
!     APW and local-orbital variables     !
!-----------------------------------------!
! energy step used for numerical calculation of energy derivatives
REAL(8) :: deapwlo
! maximum allowable APW order
INTEGER, parameter :: maxapword=4
! APW order
INTEGER :: apword(0:maxlapw,maxspecies)
! maximum of apword over all angular momenta and species
INTEGER :: apwordmax
! total number of APW coefficients (l, m and order) for each species
INTEGER :: lmoapw(maxspecies)
! polynomial order used for APW radial derivatives
INTEGER :: npapw
! APW initial linearisation energies
REAL(8) :: apwe0(maxapword,0:maxlapw,maxspecies)
!
! APW linearisation energies
REAL(8), ALLOCATABLE :: apwe(:,:,:)
! APW derivative order
INTEGER :: apwdm(maxapword,0:maxlapw,maxspecies)
! apwve is .true. if the linearisation energies are allowed to vary
LOGICAL :: apwve(maxapword,0:maxlapw,maxspecies)
! APW radial functions
REAL(8), ALLOCATABLE :: apwfr(:,:,:,:,:)
! derivate of radial functions at the muffin-tin surface
REAL(8), ALLOCATABLE :: apwdfr(:,:,:)
! maximum number of local-orbitals
INTEGER, parameter :: maxlorb=200
! maximum allowable local-orbital order
INTEGER, parameter :: maxlorbord=5
! number of local-orbitals
INTEGER :: nlorb(maxspecies)
! maximum nlorb over all species
INTEGER :: nlomax
! total number of local-orbitals
INTEGER :: nlotot
! local-orbital order
INTEGER :: lorbord(maxlorb,maxspecies)
! maximum lorbord over all species
INTEGER :: lorbordmax
! polynomial order used for local-orbital radial derivatives
INTEGER :: nplorb
! local-orbital angular momentum
INTEGER :: lorbl(maxlorb,maxspecies)
! maximum lorbl over all species
INTEGER lolmax
! (lolmax+1)^2
INTEGER lolmmax
! local-orbital initial energies
REAL(8) lorbe0(maxlorbord,maxlorb,maxspecies)
!
! local-orbital energies
REAL(8), ALLOCATABLE :: lorbe(:,:,:)
!
! local-orbital derivative order
INTEGER lorbdm(maxlorbord,maxlorb,maxspecies)
!
! lorbve is .true. if the linearisation energies are allowed to vary
LOGICAL lorbve(maxlorbord,maxlorb,maxspecies)
!
! local-orbital radial functions
REAL(8), ALLOCATABLE :: lofr(:,:,:,:)
!
! band energy search tolerance
REAL(8) epsband
! maximum allowed change in energy during band energy search; enforced only if
! default energy is less than zero
REAL(8) demaxbnd
! minimum default linearisation energy over all APWs and local-orbitals
REAL(8) e0min
! if autolinengy is .true. THEN  the fixed linearisation energies are set to the
! Fermi energy minus dlefe
LOGICAL autolinengy
! difference between linearisation and Fermi energies when autolinengy is .true.
REAL(8) dlefe
! lorbcnd is .true. if conduction state local-orbitals should be added
LOGICAL lorbcnd
! conduction state local-orbital order
INTEGER lorbordc
! excess order of the APW and local-orbital functions
INTEGER nxoapwlo
! excess local orbitals
INTEGER nxlo

end module
