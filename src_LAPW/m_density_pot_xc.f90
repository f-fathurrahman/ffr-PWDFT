module m_density_pot_xc

!---------------------------------------------------------------!
!     density, potential and exchange-correlation variables     !
!---------------------------------------------------------------!
! exchange-correlation functional type
integer xctype(3)
! exchange-correlation functional description
character(512) xcdescr
! exchange-correlation functional spin requirement
integer xcspin
! exchange-correlation functional density gradient requirement
integer xcgrad
! small constant used to stabilise non-collinear GGA
real(8) dncgga
! muffin-tin and interstitial charge density
real(8), allocatable :: rhomt(:,:),rhoir(:)
! trhonorm is .true. if the density is to be normalised after every iteration
logical trhonorm
! muffin-tin and interstitial magnetisation vector field
real(8), allocatable :: magmt(:,:,:),magir(:,:)
! tcden is .true. if the current density is to be calculated
logical tcden
! muffin-tin and interstitial current density vector field
real(8), allocatable :: cdmt(:,:,:),cdir(:,:)
! amount of smoothing to be applied to the exchange-correlation potentials and
! magnetic field
integer msmooth
! muffin-tin and interstitial Coulomb potential
real(8), allocatable :: vclmt(:,:),vclir(:)
! Poisson solver pseudocharge density constant
integer npsd
! lmaxo+npsd+1
integer lnpsd
! muffin-tin and interstitial exchange energy density
real(8), allocatable :: exmt(:,:),exir(:)
! muffin-tin and interstitial correlation energy density
real(8), allocatable :: ecmt(:,:),ecir(:)
! muffin-tin and interstitial exchange-correlation potential
real(8), allocatable :: vxcmt(:,:),vxcir(:)
! constant part of exchange-correlation potential
real(8) vxc0
! muffin-tin and interstitial Kohn-Sham effective potential
real(8), allocatable :: vsmt(:,:),vsir(:)
! G-space interstitial Kohn-Sham effective potential
complex(8), allocatable :: vsig(:)
! muffin-tin and interstitial exchange-correlation magnetic field
real(8), allocatable :: bxcmt(:,:,:),bxcir(:,:)
! muffin-tin and interstitial magnetic dipole field
real(8), allocatable :: bdmt(:,:,:),bdir(:,:)
! tbdip is .true. if the spin and current dipole fields are to be added to the
! Kohn-Sham magnetic field
logical tbdip
! muffin-tin Kohn-Sham effective magnetic field in spherical coordinates and on
! a coarse radial mesh
real(8), allocatable :: bsmt(:,:,:)
! interstitial Kohn-Sham effective magnetic field
real(8), allocatable :: bsir(:,:)
! nosource is .true. if the field is to be made source-free
logical nosource
! tssxc is .true. if scaled spin exchange-correlation (SSXC) is to be used
logical tssxc
! SSXC scaling factor
real(8) ssxc
! spin-orbit coupling radial function
real(8), allocatable :: socfr(:,:)
! kinetic energy density
real(8), allocatable :: taumt(:,:,:),tauir(:,:)
! core kinetic energy density
real(8), allocatable :: taucr(:,:,:)
! taudft is .true. if meta-GGA is to be treated as a tau-DFT functional
logical taudft
! tau-DFT exchange-correlation potential
real(8), allocatable :: wxcmt(:,:),wxcir(:)
! tau-DFT Kohn-Sham potential
real(8), allocatable :: wsmt(:,:),wsir(:)
! Tran-Blaha '09 constant c [Phys. Rev. Lett. 102, 226401 (2009)]
real(8) c_tb09
! tc_tb09 is .true. if the Tran-Blaha constant has been read in
logical tc_tb09
! if trdstate is .true. the density and potential can be read from STATE.OUT
logical trdstate
! temperature in degrees Kelvin
real(8) tempk
! maximum number of iterations used for inverting the Kohn-Sham equations
integer maxitksi
! step size used for inverting the Kohn-Sham equations
real(8) tauksi

end module
