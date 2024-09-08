
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: init0
! !INTERFACE:
subroutine init0
! !USES:
use modmain
use modxcifc
use moddftu
use modtddft
use modphonon
use modulr
use modtest
use modvars
use modmpi
use modomp
! !DESCRIPTION:
!   Performs basic consistency checks as well as allocating and initialising
!   global variables not dependent on the $k$-point set.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ist
integer nr,l,m,lm,i
real(8) rsum,t1
real(8) ts0,ts1

!-------------------------------!
!     zero timing variables     !
!-------------------------------!
timeinit=0.d0
timemat=0.d0
timefv=0.d0
timesv=0.d0
timerho=0.d0
timepot=0.d0
timefor=0.d0
call timesec(ts0)

!------------------------------------!
!     angular momentum variables     !
!------------------------------------!
if (lmaxo.gt.lmaxapw) then
  write(*,*)
  write(*,'("Error(init0): lmaxo > lmaxapw : ",2I8)') lmaxo,lmaxapw
  write(*,*)
  stop
end if
lmaxi=min(lmaxi,lmaxo)
lmmaxapw=(lmaxapw+1)**2
lmmaxi=(lmaxi+1)**2
lmmaxo=(lmaxo+1)**2
! check DOS lmax is within range
lmaxdos=min(lmaxdos,lmaxo)
! index to (l,m) pairs
if (allocated(idxlm)) deallocate(idxlm)
allocate(idxlm(0:lmaxapw,-lmaxapw:lmaxapw))
if (allocated(idxil)) deallocate(idxil)
allocate(idxil(lmmaxapw))
if (allocated(idxim)) deallocate(idxim)
allocate(idxim(lmmaxapw))
lm=0
do l=0,lmaxapw
  do m=-l,l
    lm=lm+1
    idxlm(l,m)=lm
    idxil(lm)=l
    idxim(lm)=m
  end do
end do
! array of i^l and (-i)^l values
if (allocated(zil)) deallocate(zil)
if (allocated(zilc)) deallocate(zilc)
allocate(zil(0:lmaxapw),zilc(0:lmaxapw))
do l=0,lmaxapw
  zil(l)=zi**l
  zilc(l)=conjg(zil(l))
end do
! write to VARIABLES.OUT
call writevars('lmaxapw',iv=lmaxapw)
call writevars('lmaxi',iv=lmaxi)
call writevars('lmaxo',iv=lmaxo)

!------------------------------------!
!     index to atoms and species     !
!------------------------------------!
natmmax=0
ias=0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=ias+1
    idxas(ia,is)=ias
    idxis(ias)=is
    idxia(ias)=ia
  end do
! maximum number of atoms over all species
  natmmax=max(natmmax,natoms(is))
end do
! total number of atoms
natmtot=ias
! number of phonon branches
nbph=3*natmtot
! write to VARIABLES.OUT
call writevars('nspecies',iv=nspecies)
call writevars('natoms',nv=nspecies,iva=natoms)
call writevars('spsymb',nv=nspecies,sva=spsymb)
call writevars('spname',nv=nspecies,sva=spname)
call writevars('spzn',nv=nspecies,rva=spzn)

!------------------------!
!     spin variables     !
!------------------------!
if (spinsprl) then
  spinpol=.true.
  spinorb=.false.
  if (any(task.eq.[51,52,53,61,62,63,700,701])) then
    write(*,*)
    write(*,'("Error(init0): spin-spirals do not work with task ",I4)') task
    write(*,*)
    stop
  end if
  if (xctype(1).lt.0) then
    write(*,*)
    write(*,'("Error(init0): spin-spirals do not work with the OEP method")')
    write(*,*)
    stop
  end if
  if (tefvit) then
    write(*,*)
    write(*,'("Error(init0): spin-spirals do not work with iterative &
     &diagonalisation")')
    write(*,*)
    stop
  end if
end if
! de-phasing required only for spin-spirals
if (.not.spinsprl) ssdph=.false.
! spin-orbit coupling, fixed spin moment, spin spirals or spin-polarised cores
! requires a spin-polarised calculation
if ((spinorb).or.(fsmtype.ne.0).or.(spinsprl).or.(spincore)) spinpol=.true.
! number of spinor components and maximum allowed occupancy
if (spinpol) then
  nspinor=2
  occmax=1.d0
else
  nspinor=1
  occmax=2.d0
end if
! number of spin-dependent first-variational functions per state and map from
! second- to first-variational spin index
if (spinsprl) then
  nspnfv=2
  jspnfv(1)=1
  jspnfv(2)=2
else
  nspnfv=1
  jspnfv(1)=1
  jspnfv(2)=1
end if
! no calculation of second-variational eigenvectors by default
tevecsv=.false.
! spin-polarised calculations require second-variational eigenvectors
if (spinpol) tevecsv=.true.
! Hartree-Fock/RDMFT/TDDFT/GW requires second-variational eigenvectors
if (any(task.eq.[5,10,170,300,460,461,600,620,630])) then
  tevecsv=.true.
end if
! get exchange-correlation functional data
call getxcdata(xctype,xcdescr,xcspin,xcgrad,hybrid,hybridc)
if ((spinpol).and.(xcspin.eq.0)) then
  write(*,*)
  write(*,'("Error(init0): requested spin-polarised run with &
   &spin-unpolarised")')
  write(*,'(" exchange-correlation functional")')
  write(*,*)
  stop
end if
! check for collinearity in the z-direction and set the dimension of the
! magnetisation and exchange-correlation vector fields
if (spinpol) then
  ndmag=1
  if ((abs(bfieldc0(1)).gt.epslat).or.(abs(bfieldc0(2)).gt.epslat)) ndmag=3
  do is=1,nspecies
    do ia=1,natoms(is)
      if ((abs(bfcmt0(1,ia,is)).gt.epslat).or. &
          (abs(bfcmt0(2,ia,is)).gt.epslat)) ndmag=3
    end do
  end do
! spin-orbit coupling is non-collinear in general
  if (spinorb) ndmag=3
! source-free fields and spin-spirals must be non-collinear
  if ((nosource).or.(spinsprl)) then
    ndmag=3
    cmagz=.false.
  end if
! force collinear magnetism along the z-axis if required
  if (cmagz) ndmag=1
else
  ndmag=0
end if
! set the non-collinear flag
if (ndmag.eq.3) then
  ncmag=.true.
else
  ncmag=.false.
end if
! check for meta-GGA with non-collinearity
if (((xcgrad.eq.3).or.(xcgrad.eq.4)).and.ncmag) then
  write(*,*)
  write(*,'("Error(init0): meta-GGA is not valid for non-collinear magnetism")')
  write(*,*)
  stop
end if
if (ncmag.or.spinorb) then
! spins are coupled in the second-variational Hamiltonian
  spcpl=.true.
else
! spins are decoupled
  spcpl=.false.
end if
! spin-polarised cores
if (.not.spinpol) spincore=.false.
if (fsmtype.ne.0) then
! set fixed spin moment effective field to zero
  bfsmc(:)=0.d0
! set muffin-tin FSM fields to zero
  if (allocated(bfsmcmt)) deallocate(bfsmcmt)
  allocate(bfsmcmt(3,natmtot))
  bfsmcmt(:,:)=0.d0
  if (mp_mpi.and.(mixtype.ne.1)) then
    write(*,*)
    write(*,'("Info(init0): mixtype changed to 1 for FSM calculation")')
  end if
  mixtype=1
end if
! number of independent spin components of the f_xc spin tensor
if (spinpol) then
  if (ncmag) then
    nscfxc=10
  else
    nscfxc=3
  end if
else
  nscfxc=1
end if
! set the magnetic fields to the initial values
bfieldc(:)=bfieldc0(:)
bfcmt(:,:,:)=bfcmt0(:,:,:)
! if reducebf < 1 then reduce the external magnetic fields immediately for
! non-self-consistent calculations or resumptions
if (reducebf.lt.1.d0-1.d-4) then
  if (all(task.ne.[0,1,2,3,5,28,200,201,350,351,360,630])) then
    bfieldc(:)=0.d0
    bfcmt(:,:,:)=0.d0
  end if
end if
! set the fixed tensor moment spatial and spin rotation matrices equal for the
! case of spin-orbit coupling; parity for spin is ignored by rotdmat
if (spinorb) then
  do i=1,ntmfix
    rtmfix(:,:,2,i)=rtmfix(:,:,1,i)
  end do
end if
! generate the fixed tensor moment density matrices if required
call gendmftm
! write to VARIABLES.OUT
call writevars('nspinor',iv=nspinor)
call writevars('ndmag',iv=ndmag)

!----------------------------------!
!     crystal structure set up     !
!----------------------------------!
! generate the reciprocal lattice vectors and unit cell volume
call reciplat(avec,bvec,omega,omegabz)
! inverse of the lattice vector matrix
call r3minv(avec,ainv)
! inverse of the reciprocal vector matrix
call r3minv(bvec,binv)
! Cartesian coordinates of the spin-spiral vector
call r3mv(bvec,vqlss,vqcss)
do is=1,nspecies
  do ia=1,natoms(is)
! map atomic lattice coordinates to [0,1)
    call r3frac(epslat,atposl(:,ia,is))
! determine atomic Cartesian coordinates
    call r3mv(avec,atposl(:,ia,is),atposc(:,ia,is))
  end do
end do
! check for overlapping muffin-tins and adjust radii if required
call checkmt
! compute the total muffin-tin volume (M. Meinert)
omegamt=0.d0
do is=1,nspecies
  omegamt=omegamt+dble(natoms(is))*(fourpi/3.d0)*rmt(is)**3
end do
! input q-vector in Cartesian coordinates
call r3mv(bvec,vecql,vecqc)
! write to VARIABLES.OUT
call writevars('avec',nv=9,rva=avec)
call writevars('bvec',nv=9,rva=bvec)
call writevars('omega',rv=omega)
do is=1,nspecies
  call writevars('atposl',l=is,nv=3*natoms(is),rva=atposl(:,:,is))
end do

!-------------------------------!
!     vector fields E and A     !
!-------------------------------!
tefield=.false.
if (sum(abs(efieldc(:))).gt.epslat) then
! no shift of the atomic positions
  tshift=.false.
! electric field vector in lattice coordinates
  call r3mv(ainv,efieldc,efieldl)
  tefield=.true.
end if
tafield=.false.
if (sum(abs(afieldc(:))).gt.epslat) then
  tafield=.true.
! A-field in lattice coordinates
  call r3mv(ainv,afieldc,afieldl)
! vector potential added in second-variational step
  tevecsv=.true.
end if
tafieldt=.false.
if (any(task.eq.[460,461,480,481])) then
! generate the time-step grid
  call gentimes
! read time-dependent A-field from file
  call readafieldt
  tafieldt=.true.
end if

!---------------------------------!
!     crystal symmetry set up     !
!---------------------------------!
call symmetry

!-----------------------!
!     radial meshes     !
!-----------------------!
nrmtmax=1
nrcmtmax=1
do is=1,nspecies
! make the muffin-tin mesh commensurate with lradstp
  nrmt(is)=nrmt(is)-mod(nrmt(is)-1,lradstp)
  nrmtmax=max(nrmtmax,nrmt(is))
! number of coarse radial mesh points
  nrcmt(is)=(nrmt(is)-1)/lradstp+1
  nrcmtmax=max(nrcmtmax,nrcmt(is))
end do
! set up atomic and muffin-tin radial meshes
call genrmesh
! number of points in packed muffin-tins
npmtmax=1
npcmtmax=1
do is=1,nspecies
  npmti(is)=lmmaxi*nrmti(is)
  npmt(is)=npmti(is)+lmmaxo*(nrmt(is)-nrmti(is))
  npmtmax=max(npmtmax,npmt(is))
  npcmti(is)=lmmaxi*nrcmti(is)
  npcmt(is)=npcmti(is)+lmmaxo*(nrcmt(is)-nrcmti(is))
  npcmtmax=max(npcmtmax,npcmt(is))
end do

!--------------------------------------!
!     charges and number of states     !
!--------------------------------------!
chgzn=0.d0
chgcrtot=0.d0
chgval=0.d0
nstspmax=0
nstcr=0
do is=1,nspecies
! nuclear charge
  chgzn=chgzn+spzn(is)*natoms(is)
! find the maximum number of atomic states
  nstspmax=max(nstspmax,nstsp(is))
! compute the electronic charge for each species, as well as the total core and
! valence charge
  spze(is)=0.d0
  chgcr(is)=0.d0
  do ist=1,nstsp(is)
    spze(is)=spze(is)+occsp(ist,is)
    if (spcore(ist,is)) then
      chgcr(is)=chgcr(is)+occsp(ist,is)
      nstcr=nstcr+2*ksp(ist,is)*natoms(is)
    else
      chgval=chgval+occsp(ist,is)*natoms(is)
    end if
  end do
  chgcrtot=chgcrtot+chgcr(is)*natoms(is)
end do
! add excess charge
chgval=chgval+chgexs
! total charge
chgtot=chgcrtot+chgval
if (chgtot.lt.1.d-8) then
  write(*,*)
  write(*,'("Error(init0): zero total charge")')
  write(*,*)
  stop
end if
! effective Wigner radius
rwigner=(3.d0/(fourpi*(chgtot/omega)))**(1.d0/3.d0)
! write to VARIABLES.OUT
call writevars('spze',nv=nspecies,rva=spze)
call writevars('chgcr',nv=nspecies,rva=chgcr)
call writevars('chgexs',rv=chgexs)
call writevars('chgval',rv=chgtot)

!-------------------------!
!     G-vector arrays     !
!-------------------------!
! determine gkmax from rgkmax
if (nspecies.eq.0) isgkmax=-2
select case(isgkmax)
case(:-4)
! use largest muffin-tin radius
  gkmax=rgkmax/maxval(rmt(1:nspecies))
case(-3)
! use smallest muffin-tin radius
  gkmax=rgkmax/minval(rmt(1:nspecies))
case(-2)
! use the fixed value of 2.0
  gkmax=rgkmax/2.d0
case(-1)
! use average muffin-tin radius
  rsum=0.d0
  do is=1,nspecies
     rsum=rsum+dble(natoms(is))*rmt(is)
  end do
  rsum=rsum/dble(natmtot)
  gkmax=rgkmax/rsum
case(1:)
! use user-specified muffin-tin radius
  if (isgkmax.le.nspecies) then
    gkmax=rgkmax/rmt(isgkmax)
  else
    write(*,*)
    write(*,'("Error(init0): isgkmax > nspecies : ",2I8)') isgkmax,nspecies
    write(*,*)
    stop
  end if
end select
! generate the G-vectors
call gengvec
! write number of G-vectors to test file
call writetest(900,'number of G-vectors',iv=ngvec)
! Poisson solver pseudocharge density constant
if (nspecies.gt.0) then
  t1=0.25d0*gmaxvr*maxval(rmt(1:nspecies))
else
  t1=0.25d0*gmaxvr*2.d0
end if
npsd=max(nint(t1),1)
lnpsd=lmaxo+npsd+1
! generate the Coulomb Green's function in G-space = fourpi / G^2
call gengclg
! compute the spherical Bessel functions j_l(|G|R_mt)
if (allocated(jlgrmt)) deallocate(jlgrmt)
allocate(jlgrmt(0:lnpsd,ngvec,nspecies))
call genjlgprmt(lnpsd,ngvec,gc,ngvec,jlgrmt)
! generate the spherical harmonics of the G-vectors
call genylmg
! allocate structure factor array for G-vectors
if (allocated(sfacg)) deallocate(sfacg)
allocate(sfacg(ngvec,natmtot))
! generate structure factors for G-vectors
call gensfacgp(ngvec,vgc,ngvec,sfacg)
! generate the smooth step function form factors
if (allocated(ffacg)) deallocate(ffacg)
allocate(ffacg(ngtot,nspecies))
do is=1,nspecies
  call genffacgp(is,gc,ffacg(:,is))
end do
! generate the characteristic function
call gencfun
! write to VARIABLES.OUT
call writevars('gmaxvr',rv=gmaxvr)
call writevars('ngridg',nv=3,iva=ngridg)
call writevars('intgv',nv=6,iva=intgv)
call writevars('ngvec',iv=ngvec)
call writevars('ivg',nv=3*ngtot,iva=ivg)
call writevars('igfft',nv=ngtot,iva=igfft)

!-------------------------!
!     atoms and cores     !
!-------------------------!
! determine the nuclear Coulomb potential
if (allocated(vcln)) deallocate(vcln)
allocate(vcln(nrspmax,nspecies))
t1=1.d0/y00
do is=1,nspecies
  nr=nrsp(is)
  call potnucl(ptnucl,nr,rsp(:,is),spzn(is),vcln(:,is))
  vcln(1:nr,is)=t1*vcln(1:nr,is)
end do
! solve the Kohn-Sham-Dirac equations for all atoms
call allatoms
! allocate core state occupancy and eigenvalue arrays and set to default
if (allocated(occcr)) deallocate(occcr)
allocate(occcr(nstspmax,natmtot))
if (allocated(evalcr)) deallocate(evalcr)
allocate(evalcr(nstspmax,natmtot))
do ias=1,natmtot
  is=idxis(ias)
  do ist=1,nstsp(is)
    occcr(ist,ias)=occsp(ist,is)
    evalcr(ist,ias)=evalsp(ist,is)
  end do
end do
! allocate core state radial wavefunction array
if (allocated(rwfcr)) deallocate(rwfcr)
allocate(rwfcr(nrspmax,2,nstspmax,natmtot))
! number of core spin channels
if (spincore) then
  nspncr=2
else
  nspncr=1
end if
! allocate core state charge density array
if (allocated(rhocr)) deallocate(rhocr)
allocate(rhocr(nrmtmax,natmtot,nspncr))

!-------------------------------------------------------------!
!     charge density, potentials and exchange-correlation     !
!-------------------------------------------------------------!
! allocate charge density arrays
if (allocated(rhomt)) deallocate(rhomt)
allocate(rhomt(npmtmax,natmtot))
if (allocated(rhoir)) deallocate(rhoir)
allocate(rhoir(ngtot))
! allocate magnetisation arrays
if (allocated(magmt)) deallocate(magmt)
if (allocated(magir)) deallocate(magir)
if (spinpol) then
  allocate(magmt(npmtmax,natmtot,ndmag))
  allocate(magir(ngtot,ndmag))
end if
! check if the current density should be calculated
if (tafield.or.(any(task.eq.[371,372,373,460,461]))) then
  tcden=.true.
end if
! allocate current density arrays
if (allocated(cdmt)) deallocate(cdmt)
if (allocated(cdir)) deallocate(cdir)
if (tcden) then
  allocate(cdmt(npmtmax,natmtot,3),cdir(ngtot,3))
end if
! Coulomb potential
if (allocated(vclmt)) deallocate(vclmt)
allocate(vclmt(npmtmax,natmtot))
if (allocated(vclir)) deallocate(vclir)
allocate(vclir(ngtot))
! exchange energy density
if (allocated(exmt)) deallocate(exmt)
allocate(exmt(npmtmax,natmtot))
if (allocated(exir)) deallocate(exir)
allocate(exir(ngtot))
! correlation energy density
if (allocated(ecmt)) deallocate(ecmt)
allocate(ecmt(npmtmax,natmtot))
if (allocated(ecir)) deallocate(ecir)
allocate(ecir(ngtot))
! exchange-correlation potential
if (allocated(vxcmt)) deallocate(vxcmt)
allocate(vxcmt(npmtmax,natmtot))
if (allocated(vxcir)) deallocate(vxcir)
allocate(vxcir(ngtot))
! effective Kohn-Sham potential
if (allocated(vsmt)) deallocate(vsmt)
allocate(vsmt(npmtmax,natmtot))
if (allocated(vsir)) deallocate(vsir)
allocate(vsir(ngtot))
if (allocated(vsig)) deallocate(vsig)
allocate(vsig(ngvec))
! exchange-correlation, dipole and Kohn-Sham effective magnetic fields
if (allocated(bxcmt)) deallocate(bxcmt)
if (allocated(bxcir)) deallocate(bxcir)
if (allocated(bdmt)) deallocate(bdmt)
if (allocated(bdir)) deallocate(bdir)
if (allocated(bsmt)) deallocate(bsmt)
if (allocated(bsir)) deallocate(bsir)
if (spinpol) then
  allocate(bxcmt(npmtmax,natmtot,ndmag),bxcir(ngtot,ndmag))
  if (tbdip) allocate(bdmt(npmtmax,natmtot,ndmag),bdir(ngtot,ndmag))
  allocate(bsmt(npcmtmax,natmtot,ndmag),bsir(ngtot,ndmag))
end if
! kinetic energy density
if (allocated(taumt)) deallocate(taumt)
if (allocated(tauir)) deallocate(tauir)
if (allocated(taucr)) deallocate(taucr)
if ((xcgrad.eq.3).or.(xcgrad.eq.4)) then
  allocate(taumt(npmtmax,natmtot,nspinor),tauir(ngtot,nspinor))
  allocate(taucr(npmtmax,natmtot,nspinor))
end if
! tau-DFT exchange-correlation and Kohn-Sham potentials
if (allocated(wxcmt)) deallocate(wxcmt)
if (allocated(wxcir)) deallocate(wxcir)
if (allocated(wsmt)) deallocate(wsmt)
if (allocated(wsir)) deallocate(wsir)
if (xcgrad.eq.4) then
  allocate(wxcmt(npmtmax,natmtot),wxcir(ngtot))
  allocate(wsmt(npcmtmax,natmtot),wsir(ngtot))
  tevecsv=.true.
end if
! spin-orbit coupling radial function
if (allocated(socfr)) deallocate(socfr)
if (spinorb) then
  allocate(socfr(nrcmtmax,natmtot))
end if
! allocate muffin-tin charge and moment arrays
if (allocated(chgcrlk)) deallocate(chgcrlk)
allocate(chgcrlk(natmtot))
if (allocated(chgmt)) deallocate(chgmt)
allocate(chgmt(natmtot))
if (allocated(mommt)) deallocate(mommt)
allocate(mommt(3,natmtot))
! check if scaled spin exchange-correlation (SSXC) should be used
if (abs(ssxc-1.d0).gt.1.d-6) then
  tssxc=.true.
else
  tssxc=.false.
end if

!-------------------------!
!     force variables     !
!-------------------------!
if (allocated(forcehf)) deallocate(forcehf)
allocate(forcehf(3,natmtot))
if (allocated(forceibs)) deallocate(forceibs)
allocate(forceibs(3,natmtot))
if (allocated(forcetot)) deallocate(forcetot)
allocate(forcetot(3,natmtot))

!-------------------------------------------------!
!     DFT+U and fixed tensor moment variables     !
!-------------------------------------------------!
if ((dftu.ne.0).or.(ftmtype.ne.0)) then
! density matrix elements in each muffin-tin
  if (allocated(dmatmt)) deallocate(dmatmt)
  allocate(dmatmt(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
! potential matrix elements in each muffin-tin
  if (allocated(vmatmt)) deallocate(vmatmt)
  allocate(vmatmt(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
! zero the potential matrix
  vmatmt(:,:,:,:,:)=0.d0
! require the potential matrix elements be calculated
  tvmatmt=.true.
! flags for non-zero muffin-tin potential matrices
  if (allocated(tvmmt)) deallocate(tvmmt)
  allocate(tvmmt(0:lmaxdm,natmtot))
  tvmmt(:,:)=.false.
! require second-variational eigenvectors
  tevecsv=.true.
end if
if (dftu.ne.0) then
! DFT+U energy for each atom
  if (allocated(engyadu)) deallocate(engyadu)
  allocate(engyadu(natmmax,ndftu))
! interpolation constants (alpha)
  if (allocated(alphadu)) deallocate(alphadu)
  allocate(alphadu(natmmax,ndftu))
! flag the muffin-tin potential matrices which are non-zero
  do i=1,ndftu
    is=idftu(1,i)
    if (is.gt.nspecies) then
      write(*,*)
      write(*,'("Error(init0): invalid species number : ",I8)') is
      write(*,*)
      stop
    end if
    l=idftu(2,i)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      tvmmt(l,ias)=.true.
    end do
  end do
end if
if (ftmtype.ne.0) then
! allocate and zero the fixed tensor moment potential array
  if (allocated(vmftm)) deallocate(vmftm)
  allocate(vmftm(lmmaxdm,nspinor,lmmaxdm,nspinor,natmtot))
  vmftm(:,:,:,:,:)=0.d0
! flag the muffin-tin potential matrices which are non-zero
  do i=1,ntmfix
    is=itmfix(1,i)
    ia=itmfix(2,i)
    ias=idxas(ia,is)
    l=itmfix(3,i)
    tvmmt(l,ias)=.true.
  end do
end if

!-----------------------!
!     miscellaneous     !
!-----------------------!
! determine nuclear radii and volumes
call nuclei
! determine the nuclear-nuclear energy
call energynn
! get smearing function description
call getsdata(stype,sdescr)
! get mixing type description
call getmixdata(mixtype,mixdescr)
! generate the spherical harmonic transform (SHT) matrices
call genshtmat
! allocate 1D plotting arrays
if (allocated(dvp1d)) deallocate(dvp1d)
allocate(dvp1d(nvp1d))
if (allocated(vplp1d)) deallocate(vplp1d)
allocate(vplp1d(3,npp1d))
if (allocated(dpp1d)) deallocate(dpp1d)
allocate(dpp1d(npp1d))
! initial self-consistent loop number
iscl=1
tlast=.false.
! set the Fermi energy to zero
efermi=0.d0
! set the temperature from the smearing width
tempk=swidth/kboltz

call timesec(ts1)
timeinit=timeinit+ts1-ts0

return
end subroutine
!EOC

