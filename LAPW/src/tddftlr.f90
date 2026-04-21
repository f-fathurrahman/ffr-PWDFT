
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tddftlr
use modmain
use modtddft
implicit none
! local variables
logical tq0
integer, parameter :: maxit=500
integer iq,ik,isym
integer nm,i,j,l,n
integer iw,it
real(8) v(3),t1,t2
complex(8) vfxcp,z1
character(256) fname
! allocatable arrays
integer(8), allocatable :: lock(:)
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: vchi0(:,:,:),vfxc(:,:,:)
complex(8), allocatable :: eps0(:,:,:),epsi(:,:,:),epsm(:,:,:)
complex(8), allocatable :: zw(:),a(:,:)
! initialise global variables
call init0
call init1
call init2
call init3
! check q-vector is commensurate with k-point grid
v(:)=dble(ngridk(:))*vecql(:)
v(:)=abs(v(:)-nint(v(:)))
if ((v(1) > epslat).or.(v(2) > epslat).or.(v(3) > epslat)) then
  write(*,*)
  write(*,'("Error(tddftlr): q-vector incommensurate with k-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" vecql : ",3G18.10)') vecql
  write(*,*)
  stop
end if
! find the equivalent reduced q-point
call findqpt(vecql,isym,iq)
! check if q=0
tq0=.false.
if (sum(abs(vecql(:))) < epslat) tq0=.true.
! read density and potentials from file
call readstate
! read Fermi energy from a file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
! generate the G+q-vectors and related quantities
allocate(vgqc(3,ngrf),gqc(ngrf),jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
call gengqrf(vecqc,vgqc,gqc,jlgqr,ylmgq,sfacgq)
deallocate(vgqc)
! generate the regularised Coulomb Green's function in G+q-space
allocate(gclgq(ngrf))
call gengclgq(.true.,iq,ngrf,gqc,gclgq)
gclgq(:)=sqrt(gclgq(:))
! matrix sizes
if (tq0) then
  nm=ngrf+2
else
  nm=ngrf
end if
! allocate local arrays
allocate(vchi0(nm,nm,nwrf),vfxc(nm,nm,nwrf))
allocate(eps0(nm,nm,nwrf),epsi(nm,nm,nwrf))

! compute v^1/2 chi0 v^1/2 (the symmetric version of v chi0)
vchi0(:,:,:)=0.d0

do ik=1,nkptnr
  write(*,'("Info(tddftlr): ",I6," of ",I6," k-points")') ik,nkptnr
! compute v^1/2 chi0 v^1/2
  call genvchi0(.true.,ik,lock,scissor,vecql,gclgq,jlgqr,ylmgq,sfacgq,nm,vchi0)
end do

! calculate symmetric epsilon = 1 - v^1/2 chi0 v^1/2
do i=1,nm
  do j=1,nm
    eps0(i,j,:)=-vchi0(i,j,:)
    epsi(i,j,:)=vchi0(i,j,:)
  end do
  eps0(i,i,:)=eps0(i,i,:)+1.d0
  epsi(i,i,:)=epsi(i,i,:)+1.d0
end do
allocate(a(nm,nm))
vfxcp=0.d0
it=0
10 continue
! compute vchi0 v^(-1/2) f_xc v^(-1/2) vchi0
call genvfxc(tq0,.true.,gclgq,nm,vchi0,eps0,epsi,vfxc)
! begin loop over frequencies
do iw=1,nwrf
! compute 1 - v^1/2 chi0 v^1/2 - v^(-1/2) f_xc v^(-1/2) vchi0
  a(:,:)=eps0(:,:,iw)-vfxc(:,:,iw)
! invert this matrix
  call zminv(nm,a)
! left multiply by v^1/2 chi0 v^1/2
  call zgemm('N','N',nm,nm,nm,zone,vchi0(:,:,iw),nm,a,nm,zzero,epsi(:,:,iw),nm)
! compute epsilon^(-1) = 1 + v^1/2 chi v^1/2
  do i=1,nm
    epsi(i,i,iw)=1.d0+epsi(i,i,iw)
  end do
end do
if (fxctype(1) == 210) then
! self-consistent bootstrap f_xc
  it=it+1
  if (it > maxit) then
    write(*,*)
    write(*,'("Error(tddftlr): bootstrap kernel failed to converge")')
    write(*,*)
    stop
  end if
  if (mod(it,10) == 0) then
    write(*,'("Info(tddftlr): done ",I4," bootstrap iterations")') it
    write(*,'(" head of matrix v.f_xc : ",2G18.10)') vfxc(1,1,1)
  end if
! check for convergence
  t1=abs(vfxcp)-abs(vfxc(1,1,1))
  vfxcp=vfxc(1,1,1)
  if (abs(t1) > 1.d-8) goto 10
else if (fxctype(1) == 211) then
! single iteration bootstrap
  it=it+1
  if (it <= 1) goto 10
end if
deallocate(gclgq,jlgqr)
deallocate(ylmgq,sfacgq)
deallocate(vchi0,vfxc)
! invert epsilon^(-1) to find epsilon and store in array eps0
do iw=1,nwrf
  eps0(:,:,iw)=epsi(:,:,iw)
  call zminv(nm,eps0(:,:,iw))
end do
! write G = G' = 0 components to file
  do l=1,noptcomp
    i=optcomp(1,l)
    j=optcomp(2,l)
    write(fname,'("EPSILON_TDDFT_",2I1,".OUT")') i,j
    open(50,file=trim(fname),form='FORMATTED')
    write(fname,'("EPSINV_TDDFT_",2I1,".OUT")') i,j
    open(51,file=trim(fname),form='FORMATTED')
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),dble(eps0(i,j,iw))
      write(51,'(2G18.10)') dble(wrf(iw)),dble(epsi(i,j,iw))
    end do
    write(50,*)
    write(51,*)
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),aimag(eps0(i,j,iw))
      write(51,'(2G18.10)') dble(wrf(iw)),aimag(epsi(i,j,iw))
    end do
    close(50)
    close(51)
  end do
! find the macroscopic part of epsilon by inverting the 3x3 head only
if (tq0) then
  allocate(epsm(3,3,nwrf))
  do iw=1,nwrf
    epsm(1:3,1:3,iw)=epsi(1:3,1:3,iw)
    call zminv(3,epsm(:,:,iw))
  end do
! write out the macroscopic components
  do l=1,noptcomp
    i=optcomp(1,l)
    j=optcomp(2,l)
    write(fname,'("EPSM_TDDFT_",2I1,".OUT")') i,j
    open(50,file=trim(fname),form='FORMATTED')
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),dble(epsm(i,j,iw))
    end do
    write(50,*)
    do iw=2,nwrf
      write(50,'(2G18.10)') dble(wrf(iw)),aimag(epsm(i,j,iw))
    end do
    close(50)
  end do
  allocate(zw(nwrf))
! output the Faraday angle components Delta delta and Delta beta
  do iw=2,nwrf
    zw(iw)=0.5d0*zi*epsm(1,2,iw)/sqrt(epsm(1,1,iw))
  end do
  open(50,file='FARADAY.OUT',form='FORMATTED')
  do iw=2,nwrf
    write(50,'(2G18.10)') dble(wrf(iw)),dble(zw(iw))
  end do
  write(50,*)
  do iw=2,nwrf
    write(50,'(2G18.10)') dble(wrf(iw)),aimag(zw(iw))
  end do
  close(50)
! output the Kerr angle
  do iw=2,nwrf
    zw(iw)=-epsm(1,2,iw)/(sqrt(epsm(1,1,iw))*(epsm(1,1,iw)-1.d0))
  end do
  open(50,file='KERR_TDDFT.OUT',form='FORMATTED')
  do iw=2,nwrf
    write(50,'(2G18.10)') dble(wrf(iw)),dble(zw(iw))*180.d0/pi
  end do
  write(50,*)
  do iw=2,nwrf
    write(50,'(2G18.10)') dble(wrf(iw)),aimag(zw(iw))*180.d0/pi
  end do
  close(50)
! output magnetic linear dichroism (MLD) spectrum
  t1=sin(thetamld)**2
  t2=sin(2.d0*thetamld)
  do iw=2,nwrf
    z1=epsm(1,1,iw)
    zw(iw)=t2*epsm(1,2,iw)/((z1-1.d0)*(z1-(t1*(z1+1.d0))))
  end do
  open(50,file='MLD.OUT',form='FORMATTED')
  do iw=2,nwrf
    write(50,'(2G18.10)') dble(wrf(iw)),dble(zw(iw))
  end do
  write(50,*)
  do iw=2,nwrf
    write(50,'(2G18.10)') dble(wrf(iw)),aimag(zw(iw))
  end do
  close(50)
  deallocate(epsm,zw)
end if

! write inverse epsilon to test file
!call writetest(320,'inverse epsilon',nv=nm*nm*nwrf,tol=1.d-2,zva=epsi)

deallocate(eps0,epsi,a)
  write(*,*)
  write(*,'("Info(tddftlr):")')
  write(*,'(" Dielectric tensor written to EPSILON_TDDFT_ij.OUT")')
  write(*,'(" Inverse written to EPSINV_TDDFT_ij.OUT")')
  if (tq0) then
    write(*,'(" Macroscopic part written to EPSM_TDDFT_ij.OUT")')
  end if
  write(*,'(" for components")')
  do l=1,noptcomp
    write(*,'("  i = ",I1,", j = ",I1)') optcomp(1:2,l)
  end do
  write(*,'(" q-vector (lattice coordinates) : ")')
  write(*,'(3G18.10)') vecql
  write(*,'(" q-vector length : ",G18.10)') gqc(1)
  if (tq0) then
    write(*,*)
    write(*,'(" Faraday angle parameters Δδ and Δβ written to FARADAY.OUT")')
    write(*,'(" MOKE Kerr angle written to KERR_TDDFT.OUT")')
    write(*,'(" Magnetic linear dichroism (MLD) spectrum written to MLD.OUT")')
  end if

deallocate(gqc)
return
end subroutine

