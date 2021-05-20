MODULE m_density_pot_xc

!---------------------------------------------------------------!
!     density, potential and exchange-correlation variables     !
!---------------------------------------------------------------!
! exchange-correlation functional type
INTEGER xctype(3)
! exchange-correlation functional description
character(512) xcdescr
! exchange-correlation functional spin requirement
INTEGER xcspin
! exchange-correlation functional density gradient requirement
INTEGER xcgrad
! small constant used to stabilise non-collinear GGA
REAL(8) dncgga
! muffin-tin and interstitial charge density
REAL(8), ALLOCATABLE :: rhomt(:,:),rhoir(:)
! trhonorm is .true. if the density is to be normalised after every iteration
LOGICAL trhonorm
! muffin-tin and interstitial magnetisation vector field
REAL(8), ALLOCATABLE :: magmt(:,:,:),magir(:,:)
! tcden is .true. if the current density is to be calculated
LOGICAL tcden
! muffin-tin and interstitial current density vector field
REAL(8), ALLOCATABLE :: cdmt(:,:,:),cdir(:,:)
! amount of smoothing to be applied to the exchange-correlation potentials and
! magnetic field
INTEGER msmooth
! muffin-tin and interstitial Coulomb potential
REAL(8), ALLOCATABLE :: vclmt(:,:),vclir(:)
! Poisson solver pseudocharge density constant
INTEGER npsd
! lmaxo+npsd+1
INTEGER lnpsd
! muffin-tin and interstitial exchange energy density
REAL(8), ALLOCATABLE :: exmt(:,:),exir(:)
! muffin-tin and interstitial correlation energy density
REAL(8), ALLOCATABLE :: ecmt(:,:),ecir(:)
! muffin-tin and interstitial exchange-correlation potential
REAL(8), ALLOCATABLE :: vxcmt(:,:),vxcir(:)
! constant part of exchange-correlation potential
REAL(8) vxc0
! muffin-tin and interstitial Kohn-Sham effective potential
REAL(8), ALLOCATABLE :: vsmt(:,:),vsir(:)
! G-space interstitial Kohn-Sham effective potential
COMPLEX(8), ALLOCATABLE :: vsig(:)
! muffin-tin and interstitial exchange-correlation magnetic field
REAL(8), ALLOCATABLE :: bxcmt(:,:,:),bxcir(:,:)
! muffin-tin and interstitial magnetic dipole field
REAL(8), ALLOCATABLE :: bdmt(:,:,:),bdir(:,:)
! tbdip is .true. if the spin and current dipole fields are to be added to the
! Kohn-Sham magnetic field
LOGICAL tbdip
! muffin-tin Kohn-Sham effective magnetic field in spherical coordinates and on
! a coarse radial mesh
REAL(8), ALLOCATABLE :: bsmt(:,:,:)
! interstitial Kohn-Sham effective magnetic field
REAL(8), ALLOCATABLE :: bsir(:,:)
! nosource is .true. if the field is to be made source-free
LOGICAL nosource
! tssxc is .true. if scaled spin exchange-correlation (SSXC) is to be used
LOGICAL tssxc
! SSXC scaling factor
REAL(8) ssxc
! spin-orbit coupling radial function
REAL(8), ALLOCATABLE :: socfr(:,:)
! kinetic energy density
REAL(8), ALLOCATABLE :: taumt(:,:,:),tauir(:,:)
! core kinetic energy density
REAL(8), ALLOCATABLE :: taucr(:,:,:)
! taudft is .true. if meta-GGA is to be treated as a tau-DFT functional
LOGICAL taudft
! tau-DFT exchange-correlation potential
REAL(8), ALLOCATABLE :: wxcmt(:,:),wxcir(:)
! tau-DFT Kohn-Sham potential
REAL(8), ALLOCATABLE :: wsmt(:,:),wsir(:)
! Tran-Blaha '09 constant c [Phys. Rev. Lett. 102, 226401 (2009)]
REAL(8) c_tb09
! tc_tb09 is .true. if the Tran-Blaha constant has been read in
LOGICAL tc_tb09
! if trdstate is .true. the density and potential can be read from STATE.OUT
LOGICAL trdstate
! temperature in degrees Kelvin
REAL(8) tempk
! maximum number of iterations used for inverting the Kohn-Sham equations
INTEGER maxitksi
! step size used for inverting the Kohn-Sham equations
REAL(8) tauksi

END MODULE 
