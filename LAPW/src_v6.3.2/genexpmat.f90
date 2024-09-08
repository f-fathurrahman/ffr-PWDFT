
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genexpmat(vpl,expmt,emat)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(in) :: expmt(npcmtmax,natmtot)
complex(8), intent(out) :: emat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,i,j,k,l
integer is,ia,ias,nrc,nrci
integer npc,ngp,ngpq,igp,ifg
real(8) vpc(3),vpql(3),vpqc(3),t1
complex(8) z1
! allocatable arrays
integer, allocatable :: igpig(:),igpqig(:)
real(8), allocatable :: vgpl(:,:),vgpc(:,:),gpc(:)
real(8), allocatable :: vgpql(:,:),vgpqc(:,:),gpqc(:)
complex(8), allocatable :: sfacgp(:,:),sfacgpq(:,:)
complex(8), allocatable :: apwalm1(:,:,:,:),apwalm2(:,:,:,:)
complex(8), allocatable :: evecfv1(:,:),evecfv2(:,:)
complex(8), allocatable :: evecsv1(:,:),evecsv2(:,:)
complex(8), allocatable :: wfmt1(:),wfmt2(:,:)
complex(8), allocatable :: zfir(:),x(:),em(:,:)
! external functions
complex(8) zfmtinp,zdotc
external zfmtinp,zdotc
! check if q-vector is zero
t1=abs(vecql(1))+abs(vecql(2))+abs(vecql(3))
if (t1.lt.epslat) then
  emat(:,:)=0.d0
  do i=1,nstsv
    emat(i,i)=1.d0
  end do
  return
end if
! allocate local arrays
allocate(igpig(ngkmax),igpqig(ngkmax))
allocate(vgpl(3,ngkmax),vgpc(3,ngkmax),gpc(ngkmax))
allocate(vgpql(3,ngkmax),vgpqc(3,ngkmax),gpqc(ngkmax))
allocate(sfacgp(ngkmax,natmtot),sfacgpq(ngkmax,natmtot))
allocate(apwalm1(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(apwalm2(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv1(nmatmax,nstfv),evecfv2(nmatmax,nstfv))
if (tevecsv) then
  allocate(evecsv1(nstsv,nstsv),evecsv2(nstsv,nstsv))
end if
allocate(wfmt1(npcmtmax),wfmt2(npcmtmax,nstfv))
allocate(zfir(ngtot),x(ngkmax),em(nstfv,nstfv))
! p-vector in Cartesian coordinates
call r3mv(bvec,vpl,vpc)
! generate the G+p-vectors
call gengkvec(ngvec,ivg,vgc,vpl,vpc,gkmax,ngkmax,ngp,igpig,vgpl,vgpc,gpc)
! generate the structure factors
call gensfacgp(ngp,vgpc,ngkmax,sfacgp)
! find the matching coefficients for k-point p
call match(ngp,vgpc,gpc,sfacgp,apwalm1)
! get the eigenvectors for k-point p
call getevecfv(filext,0,vpl,vgpl,evecfv1)
! p+q-vector in lattice coordinates
vpql(:)=vpl(:)+vecql(:)
! p+q-vector in Cartesian coordinates
call r3mv(bvec,vpql,vpqc)
! generate the G+p+q-vectors
call gengkvec(ngvec,ivg,vgc,vpql,vpqc,gkmax,ngkmax,ngpq,igpqig,vgpql,vgpqc,gpqc)
! generate the structure factors
call gensfacgp(ngpq,vgpqc,ngkmax,sfacgpq)
! find the matching coefficients for k-point p+q
call match(ngpq,vgpqc,gpqc,sfacgpq,apwalm2)
! get the eigenvectors for k-point p+q
call getevecfv(filext,0,vpql,vgpql,evecfv2)
! set the first-variational matrix element array to zero
em(:,:)=0.d0
!------------------------------------!
!     muffin-tin matrix elements     !
!------------------------------------!
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist=1,nstfv
! calculate the wavefunction for k-point p+q
      call wavefmt(lradstp,ias,ngpq,apwalm2(:,:,:,ias),evecfv2(:,ist),wfmt1)
! convert from spherical harmonics to spherical coordinates
      call zbsht(nrc,nrci,wfmt1,wfmt2(:,ist))
! multiply by exp(-iq.r) (conjugate because zfmtinp conjugates first function)
      wfmt1(1:npc)=conjg(expmt(1:npc,ias))*wfmt2(1:npc,ist)
! convert from spherical coordinates to spherical harmonics
      call zfsht(nrc,nrci,wfmt1,wfmt2(:,ist))
    end do
    do jst=1,nstfv
! calculate the wavefunction for k-point p
      call wavefmt(lradstp,ias,ngp,apwalm1(:,:,:,ias),evecfv1(:,jst),wfmt1)
      do ist=1,nstfv
        em(ist,jst)=em(ist,jst)+zfmtinp(nrc,nrci,wrcmt(:,is),wfmt2(:,ist),wfmt1)
      end do
    end do
! end loops over atoms and species
  end do
end do
!--------------------------------------!
!     interstitial matrix elements     !
!--------------------------------------!
! compute interstitial wavefunctions for k-point p
do jst=1,nstfv
  zfir(:)=0.d0
  do igp=1,ngp
    ifg=igfft(igpig(igp))
    zfir(ifg)=evecfv1(igp,jst)
  end do
! Fourier transform wavefunction to real-space
  call zfftifc(3,ngridg,1,zfir)
! multiply with the characteristic function
  zfir(:)=zfir(:)*cfunir(:)
! Fourier transform back to G-space
  call zfftifc(3,ngridg,-1,zfir)
! store as wavefunction with G+p+q index
  do igp=1,ngpq
    ifg=igfft(igpqig(igp))
    x(igp)=zfir(ifg)
  end do
! add to the first-variational matrix elements
  do ist=1,nstfv
    em(ist,jst)=em(ist,jst)+zdotc(ngpq,evecfv2(:,ist),1,x,1)
  end do
end do
!-------------------------------------------!
!     second-variational matrix elements    !
!-------------------------------------------!
if (tevecsv) then
! get the second-variational eigenvectors
  call getevecsv(filext,0,vpl,evecsv1)
  call getevecsv(filext,0,vpql,evecsv2)
  do i=1,nstsv
    do j=1,nstsv
      z1=0.d0
      k=0
      do ispn=1,nspinor
        do ist=1,nstfv
          k=k+1
          l=(ispn-1)*nstfv
          do jst=1,nstfv
            l=l+1
            z1=z1+em(ist,jst)*conjg(evecsv2(k,i))*evecsv1(l,j)
          end do
        end do
      end do
      emat(i,j)=z1
    end do
  end do
else
  emat(:,:)=em(:,:)
end if
deallocate(igpig,igpqig,vgpl,vgpc,gpc)
deallocate(vgpql,vgpqc,gpqc,sfacgp,sfacgpq)
deallocate(apwalm1,apwalm2,evecfv1,evecfv2)
if (tevecsv) deallocate(evecsv1,evecsv2)
deallocate(wfmt1,wfmt2,zfir,x,em)
return
end subroutine

