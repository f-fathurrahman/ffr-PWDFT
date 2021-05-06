SUBROUTINE read_input()
  USE m_atoms, ONLY: natoms, xctsp, atposl, atposc, ecvcut, esccut, molecule, nspecies, &
               primcell, ptnucl, rndatposc, sppath
  !USE m_apwlo
  USE m_muffin_tins
  USE m_lattice
  USE m_dos_optics_response
  USE m_plotting
  USE m_mixing
  USE m_density_pot_xc
  USE m_misc
  USE m_symmetry
  USE m_force_stress
  USE m_hamiltonian
  USE m_states
  USE m_spin
  USE m_core_states
  USE m_gvectors
  USE m_gkvectors
  USE m_qpoints
  USE m_kpoints
  USE m_convergence
  USE m_oep_hf
  USE m_charge_moment_current
  USE m_electric_vector_pot, ONLY: afieldc, efieldc
  IMPLICIT NONE 
  ! local variables
  logical :: lv
  INTEGER :: is,ia,ias,ios
  INTEGER :: i,j,k,l,n,p
  REAL(8) :: sc,sc1,sc2,sc3
  REAL(8) :: scu,scu1,scu2,scu3
  REAL(8) :: solscf,zn,a,b
  REAL(8) :: axang(4),rot(3,3)
  REAL(8) :: v1(3),v2(3),t1
  character(256) :: block,symb,str

  !------------------------!
  !     default values     !
  !------------------------!
  ntasks=0
  avec(:,:)=0.d0
  avec(1,1)=1.d0
  avec(2,2)=1.d0
  avec(3,3)=1.d0
  sc=1.d0
  sc1=1.d0
  sc2=1.d0
  sc3=1.d0
  epslat=1.d-6
  primcell=.false.
  tshift=.true.
  ngridk(:)=1
  vkloff(:)=0.d0
  autokpt=.false.
  radkpt=40.d0
  reducek=1
  ngridq(:)=-1
  reduceq=1
  rgkmax=7.d0
  gmaxvr=12.d0
  lmaxapw=8
  lmaxo=6
  lmaxi=1
  fracinr=0.01d0
  trhonorm=.true.
  xctype(1)=3
  xctype(2:3)=0
  xctsp(1)=3
  xctsp(2:3)=0
  stype=3
  swidth=0.001d0
  autoswidth=.false.
  mstar=10.d0
  epsocc=1.d-8
  epschg=1.d-3
  nempty0=4.d0
  maxscl=200
  mixtype=3
  amixpm(1)=0.05d0
  amixpm(2)=1.d0
  ! Broyden parameters recommended by M. Meinert
  mixsdb=5
  broydpm(1)=0.4d0
  broydpm(2)=0.15d0
  epspot=1.d-6
  epsengy=1.d-4
  epsforce=5.d-3
  epsstress=1.d-3
  molecule=.false.
  nspecies=0
  natoms(:)=0
  atposl(:,:,:)=0.d0
  atposc(:,:,:)=0.d0
  bfcmt0(:,:,:)=0.d0
  sppath=''
  scrpath=''
  nvp1d=2
  IF(allocated(vvlp1d)) DEALLOCATE(vvlp1d)
  ALLOCATE(vvlp1d(3,nvp1d))
  vvlp1d(:,1)=0.d0
  vvlp1d(:,2)=1.d0
  npp1d=200
  vclp2d(:,:)=0.d0
  vclp2d(1,1)=1.d0
  vclp2d(2,2)=1.d0
  np2d(:)=40
  vclp3d(:,:)=0.d0
  vclp3d(1,1)=1.d0
  vclp3d(2,2)=1.d0
  vclp3d(3,3)=1.d0
  np3d(:)=20
  nwplot=500
  ngrkf=100
  nswplot=1
  wplot(1)=-0.5d0
  wplot(2)=0.5d0
  dosocc=.false.
  dosmsum=.false.
  dosssum=.false.
  lmirep=.true.
  spinpol=.false.
  spinorb=.false.
  socscf=1.d0
  atpopt=1
  maxatpstp=200
  tau0atp=0.25d0
  deltast=0.001d0
  latvopt=0
  maxlatvstp=30
  tau0latv=0.2d0
  lradstp=4
  chgexs=0.d0
  scissor=0.d0
  noptcomp=1
  optcomp(:,1)=1
  intraband=.false.
  evaltol=-1.d0
  epsband=1.d-12
  demaxbnd=2.5d0
  autolinengy=.false.
  dlefe=-0.1d0
  deapwlo=0.05d0
  bfieldc0(:)=0.d0
  efieldc(:)=0.d0
  afieldc(:)=0.d0
  fsmtype=0
  momfix(:)=0.d0
  mommtfix(:,:,:)=1.d6

  rmtdelta=0.05d0
  isgkmax=-1
  symtype=1
  notelns=0
  tforce=.false.
  maxitoep=200

  nkstlist=1
  kstlist(:,1)=1
  vklem(:)=0.d0
  deltaem=0.025d0
  ndspem=1
  nosource=.false.
  spinsprl=.false.
  ssdph=.true.
  vqlss(:)=0.d0
  nwrite=0

  reducebf=1.d0
  ptnucl=.true.
  tefvr=.true.
  tefvit=.false.
  minitefv=6
  maxitefv=4
  befvit=0.25d0
  epsefvit=1.d-5
  vecql(:)=0.d0
  sqados(1:2)=0.d0
  sqados(3)=1.d0

  spincore=.false.
  solscf=1.d0
  emaxelnes=-1.2d0
  wsfac(1)=-1.d6; wsfac(2)=1.d6
  
  hybrid=.false.
  hybridc=1.d0
  ecvcut=-3.5d0
  esccut=-0.4d0
  
  rndatposc=0.d0
  rndbfcmt=0.d0
  rndavec=0.d0
  c_tb09=0.d0
  tc_tb09=.false.
  lorbcnd=.false.
  lorbordc=3
  nrmtscf=1.d0
  lmaxdos=3
  msmooth=0
  npmae0=-1
  wrtvars=.false.
  cmagz=.false.
  axang(:)=0.d0
  dncgga=1.d-8

  nxoapwlo=0
  nxlo=0

  evtype=1

  maxitksi=200
  tauksi=0.002d0
  tfav0=.true.
  
  rmtall=-1.d0
  taudft=.false.

  ! read in atomic species data
  CALL read_species()
  
  RETURN 

!CONTAINS 
!
!  SUBROUTINE addstr(slist)
!    IMPLICIT NONE 
!    ! arguments
!    character(256), intent(inout), ALLOCATABLE :: slist(:)
!    ! ALLOCATABLE arrays
!    character(256), ALLOCATABLE :: stmp(:)
!    n=size(slist)
!    ALLOCATE(stmp(n))
!    stmp(1:n)=slist(1:n)
!    DEALLOCATE(slist)
!    ALLOCATE(slist(n+1))
!    slist(1:n)=stmp(1:n)
!    slist(n+1)=str
!    DEALLOCATE(stmp)
!    RETURN 
!  END SUBROUTINE 

END SUBROUTINE 