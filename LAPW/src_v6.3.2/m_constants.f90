module m_constants

!-----------------------------!
!     numerical constants     !
!-----------------------------!
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: twopi=6.2831853071795864769d0
real(8), parameter :: fourpi=12.566370614359172954d0
! spherical harmonic for l=m=0
real(8), parameter :: y00=0.28209479177387814347d0
! complex constants
complex(8), parameter :: zzero=(0.d0,0.d0)
complex(8), parameter :: zone=(1.d0,0.d0)
complex(8), parameter :: zi=(0.d0,1.d0)
! array of i^l and (-i)^l values
complex(8), allocatable :: zil(:),zilc(:)
! Pauli spin matrices:
! sigma_x = ( 0  1 )   sigma_y = ( 0 -i )   sigma_z = ( 1  0 )
!           ( 1  0 )             ( i  0 )             ( 0 -1 )
complex(8) sigmat(2,2,3)
data sigmat / (0.d0,0.d0), (1.d0,0.d0), (1.d0,0.d0), (0.d0,0.d0), &
              (0.d0,0.d0), (0.d0,1.d0),(0.d0,-1.d0), (0.d0,0.d0), &
              (1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0),(-1.d0,0.d0) /
! Planck constant in SI units (exact, CODATA 2018)
real(8), parameter :: h_si=6.62607015d-34
! reduced Planck constant in SI units
real(8), parameter :: hbar_si=h_si/twopi
! speed of light in SI units (exact, CODATA 2018)
real(8), parameter :: sol_si=299792458d0
! speed of light in atomic units (=1/alpha) (CODATA 2018)
real(8), parameter :: sol=137.035999084d0
! scaled speed of light
real(8) solsc
! Hartree in SI units (CODATA 2018)
real(8), parameter :: ha_si=4.3597447222071d-18
! Hartree in eV (CODATA 2018)
real(8), parameter :: ha_ev=27.211386245988d0
! Hartree in inverse meters
real(8), parameter :: ha_im=ha_si/(h_si*sol_si)
! Boltzmann constant in SI units (exact, CODATA 2018)
real(8), parameter :: kb_si=1.380649d-23
! Boltzmann constant in Hartree/kelvin
real(8), parameter :: kboltz=kb_si/ha_si
! electron charge in SI units (exact, CODATA 2018)
real(8), parameter :: e_si=1.602176634d-19
! Bohr radius in SI units (CODATA 2018)
real(8), parameter :: br_si=0.529177210903d-10
! Bohr radius in Angstroms
real(8), parameter :: br_ang=br_si*1.d10
! atomic unit of magnetic flux density in SI
real(8), parameter :: b_si=hbar_si/(e_si*br_si**2)
! atomic unit of electric field in SI
real(8), parameter :: ef_si=ha_si/(e_si*br_si)
! atomic unit of time in SI
real(8), parameter :: t_si=hbar_si/ha_si
! electron g-factor (CODATA 2018)
real(8), parameter :: gfacte=2.00231930436256d0
! electron mass in SI (CODATA 2018)
real(8), parameter :: em_si=9.1093837015d-31
! atomic mass unit in SI (CODATA 2018)
real(8), parameter :: amu_si=1.66053906660d-27
! atomic mass unit in electron masses
real(8), parameter :: amu=amu_si/em_si

end module
