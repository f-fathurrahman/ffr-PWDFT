
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhoinit
! !INTERFACE:
subroutine rhoinit
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Initialises the crystal charge density. Inside the muffin-tins it is set to
!   the spherical atomic density. In the interstitial region it is taken to be
!   constant such that the total charge is correct. Requires that the atomic
!   densities have already been calculated.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer lmax,is,ia,ias,i
integer nr,nri,nro,nrs,ir
integer nrc,nrci,irco,irc
integer l,m,lm,ig,ifg,nthd
real(8) x,t1,t2
complex(8) z1,z2,z3
! allocatable arrays
real(8), allocatable :: jl(:,:),ffg(:),wr(:),fr(:)
complex(8), allocatable :: zfmt(:),zfft(:)
lmax=min(lmaxi,1)
! zero the charge density arrays
rhomt(:,:)=0.d0
rhoir(:)=0.d0
! compute the superposition of all the atomic density tails
allocate(zfft(ngtot))
zfft(:)=0.d0
call holdthd(nspecies,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ffg,wr,fr,nr,nrs,nro,ig) &
!$OMP PRIVATE(ir,t1,t2,x,ia,ias,ifg) &
!$OMP NUM_THREADS(nthd)
allocate(ffg(ngvec),wr(nrspmax),fr(nrspmax))
!$OMP DO
do is=1,nspecies
  nr=nrmt(is)
  nrs=nrsp(is)
  nro=nrs-nr+1
! determine the weights for the radial integral
  call wsplint(nro,rsp(nr,is),wr(nr))
  do ig=1,ngvec
    t1=gc(ig)
! spherical bessel function j_0(x) times the atomic density tail
    if (t1.gt.epslat) then
      t2=1.d0/t1
      do ir=nr,nrs
        x=t1*rsp(ir,is)
        fr(ir)=t2*sin(x)*rhosp(ir,is)*rsp(ir,is)
      end do
    else
      fr(nr:nrs)=rhosp(nr:nrs,is)*rsp(nr:nrs,is)**2
    end if
    t1=dot_product(wr(nr:nrs),fr(nr:nrs))
! apply low-pass filter
    t1=t1*exp(-4.d0*(gc(ig)/gmaxvr)**2)
    ffg(ig)=(fourpi/omega)*t1
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!$OMP CRITICAL(rhoinit_)
    do ig=1,ngvec
      ifg=igfft(ig)
      zfft(ifg)=zfft(ifg)+ffg(ig)*conjg(sfacg(ig,ias))
    end do
!$OMP END CRITICAL(rhoinit_)
  end do
end do
!$OMP END DO
deallocate(ffg,wr,fr)
!$OMP END PARALLEL
call freethd(nthd)
! compute the tails in each muffin-tin
call holdthd(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(jl,zfmt,is,nrc,nrci) &
!$OMP PRIVATE(irco,ig,ifg,irc,x) &
!$OMP PRIVATE(z1,z2,z3,lm,l,m,i) &
!$OMP NUM_THREADS(nthd)
allocate(jl(0:lmax,nrcmtmax),zfmt(npcmtmax))
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  irco=nrci+1
  zfmt(1:npcmt(is))=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    do irc=1,nrc
      x=gc(ig)*rcmt(irc,is)
      call sbessel(lmax,x,jl(:,irc))
    end do
    z1=fourpi*zfft(ifg)*sfacg(ig,ias)
    lm=0
    do l=0,lmax
      z2=z1*zil(l)
      do m=-l,l
        lm=lm+1
        z3=z2*conjg(ylmg(lm,ig))
        i=lm
        do irc=1,nrci
          zfmt(i)=zfmt(i)+jl(l,irc)*z3
          i=i+lmmaxi
        end do
        do irc=irco,nrc
          zfmt(i)=zfmt(i)+jl(l,irc)*z3
          i=i+lmmaxo
        end do
      end do
    end do
  end do
  call ztorfmt(nrc,nrci,zfmt,rhomt(:,ias))
end do
!$OMP END DO
deallocate(jl,zfmt)
!$OMP END PARALLEL
call freethd(nthd)
! convert the density from a coarse to a fine radial mesh
call rfmtctof(rhomt)
! add the atomic charge density and the excess charge in each muffin-tin
t1=chgexs/omega
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  i=1
  do ir=1,nri
    t2=(t1+rhosp(ir,is))/y00
    rhomt(i,ias)=rhomt(i,ias)+t2
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    t2=(t1+rhosp(ir,is))/y00
    rhomt(i,ias)=rhomt(i,ias)+t2
    i=i+lmmaxo
  end do
end do
! interstitial density determined from the atomic tails and excess charge
call zfftifc(3,ngridg,1,zfft)
do ir=1,ngtot
  rhoir(ir)=dble(zfft(ir))+t1
! make sure that the density is always positive
  if (rhoir(ir).lt.1.d-10) rhoir(ir)=1.d-10
end do
deallocate(zfft)
return
end subroutine
!EOC

