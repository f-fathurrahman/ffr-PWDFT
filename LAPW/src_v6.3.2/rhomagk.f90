
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhomagk
! !INTERFACE:
subroutine rhomagk(ngp,igpig,lock,wppt,occsvp,apwalm,evecfv,evecsv)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer(nspnfv))
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv))
!   lock   : OpenMP lock for each atom (in,integer(natmtot))
!   wppt   : weight of input p-point (in,real)
!   occsvp : occupation number for each state (in,real(nstsv))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Generates the partial valence charge density and magnetisation from the
!   eigenvectors at a particular $k$-point. In the muffin-tin region, the
!   wavefunction is obtained in terms of its $(l,m)$-components from both the
!   APW and local-orbital functions. Using a backward spherical harmonic
!   transform (SHT), the wavefunction is converted to real-space and the density
!   obtained from its modulus squared. A similar proccess is used for the
!   intersitial density in which the wavefunction in real-space is obtained from
!   a Fourier transform of the APW functions. See routines {\tt wavefmt},
!   {\tt genshtmat} and {\tt eveqn}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Removed conversion to spherical harmonics, January 2009 (JKD)
!   Partially de-phased the muffin-tin magnetisation for spin-spirals,
!    February 2009 (FC, FB & LN)
!   Optimisations, July 2010 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
integer(8), intent(in) :: lock(natmtot)
real(8), intent(in) :: wppt,occsvp(nstsv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
! local variables
integer ispn,jspn,ist
integer is,ia,ias,i,j
integer nrc,nrci,npc
integer igp,ifg,nthd
real(8) wo,t1
real(8) ts0,ts1
complex(8) zq(2),z1
! automatic arrays
logical done(nstfv,nspnfv)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:),wfmt2(:),wfmt3(:,:),wfir(:,:)
call timesec(ts0)
!----------------------------------------------!
!     muffin-tin density and magnetisation     !
!----------------------------------------------!
call holdthd(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,wfmt2,wfmt3) &
!$OMP PRIVATE(is,nrc,nrci,npc) &
!$OMP PRIVATE(ia,t1,zq,done,i,j) &
!$OMP PRIVATE(wo,ispn,jspn,ist,z1) &
!$OMP NUM_THREADS(nthd)
if (tevecsv) allocate(wfmt1(npcmtmax,nstfv,nspnfv))
allocate(wfmt2(npcmtmax),wfmt3(npcmtmax,nspinor))
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
! de-phasing factor for spin-spirals
  if (ssdph) then
    ia=idxia(ias)
    t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
    zq(1)=cmplx(cos(t1),sin(t1),8)
    zq(2)=conjg(zq(1))
  end if
  done(:,:)=.false.
  do j=1,nstsv
    wo=occsvp(j)
    if (abs(wo).lt.epsocc) cycle
    wo=wo*wppt
    if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
      i=0
      do ispn=1,nspinor
        jspn=jspnfv(ispn)
        wfmt3(1:npc,ispn)=0.d0
        do ist=1,nstfv
          i=i+1
          z1=evecsv(i,j)
          if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
            if (ssdph) z1=z1*zq(ispn)
            if (.not.done(ist,jspn)) then
              call wavefmt(lradstp,ias,ngp(jspn),apwalm(:,:,:,ias,jspn), &
               evecfv(:,ist,jspn),wfmt2)
! convert to spherical coordinates
              call zbsht(nrc,nrci,wfmt2,wfmt1(:,ist,jspn))
              done(ist,jspn)=.true.
            end if
! add to spinor wavefunction
            call zaxpy(npc,z1,wfmt1(:,ist,jspn),1,wfmt3(:,ispn),1)
          end if
        end do
      end do
    else
! spin-unpolarised wavefunction
      call wavefmt(lradstp,ias,ngp,apwalm(:,:,:,ias,1),evecfv(:,j,1),wfmt2)
! convert to spherical coordinates
      call zbsht(nrc,nrci,wfmt2,wfmt3)
    end if
! add to density and magnetisation
    call omp_set_lock(lock(ias))
    if (spinpol) then
! spin-polarised
      if (ncmag) then
! non-collinear
        call rmk1(npc,wo,wfmt3,wfmt3(:,2),rhomt(:,ias),magmt(:,ias,1), &
         magmt(:,ias,2),magmt(:,ias,3))
      else
! collinear
        call rmk2(npc,wo,wfmt3,wfmt3(:,2),rhomt(:,ias),magmt(:,ias,1))
      end if
    else
! spin-unpolarised
      call rmk3(npc,wo,wfmt3,rhomt(:,ias))
    end if
    call omp_unset_lock(lock(ias))
  end do
! end loop over atoms
end do
!$OMP END DO
if (tevecsv) deallocate(wfmt1)
deallocate(wfmt2,wfmt3)
!$OMP END PARALLEL
call freethd(nthd)
!------------------------------------------------!
!     interstitial density and magnetisation     !
!------------------------------------------------!
call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir,wo,i,ispn,jspn) &
!$OMP PRIVATE(ist,z1,igp,ifg) &
!$OMP NUM_THREADS(nthd)
allocate(wfir(ngtot,nspinor))
!$OMP DO
do j=1,nstsv
  wo=occsvp(j)
  if (abs(wo).lt.epsocc) cycle
  wo=wo*wppt/omega
  wfir(:,:)=0.d0
  if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
    i=0
    do ispn=1,nspinor
      jspn=jspnfv(ispn)
      do ist=1,nstfv
        i=i+1
        z1=evecsv(i,j)
        if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
          do igp=1,ngp(jspn)
            ifg=igfft(igpig(igp,jspn))
            wfir(ifg,ispn)=wfir(ifg,ispn)+z1*evecfv(igp,ist,jspn)
          end do
        end if
      end do
    end do
  else
! spin-unpolarised wavefunction
    do igp=1,ngp(1)
      ifg=igfft(igpig(igp,1))
      wfir(ifg,1)=evecfv(igp,j,1)
    end do
  end if
! Fourier transform wavefunction to real-space
  do ispn=1,nspinor
    call zfftifc(3,ngridg,1,wfir(:,ispn))
  end do
! add to density and magnetisation
!$OMP CRITICAL(rhomagk_)
  if (spinpol) then
! spin-polarised
    if (ncmag) then
! non-collinear
      call rmk1(ngtot,wo,wfir,wfir(:,2),rhoir,magir,magir(:,2),magir(:,3))
    else
! collinear
      call rmk2(ngtot,wo,wfir,wfir(:,2),rhoir,magir)
    end if
  else
! spin-unpolarised
    call rmk3(ngtot,wo,wfir,rhoir)
  end if
!$OMP END CRITICAL(rhomagk_)
end do
!$OMP END DO
deallocate(wfir)
!$OMP END PARALLEL
call freethd(nthd)
call timesec(ts1)
!$OMP ATOMIC
timerho=timerho+ts1-ts0
return

contains

pure subroutine rmk1(n,wo,wf1,wf2,rho,mag1,mag2,mag3)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag1(n),mag2(n),mag3(n)
! local variables
integer i
real(8) wo2,t1,t2
complex(8) z1,z2
wo2=2.d0*wo
do i=1,n
  z1=wf1(i); z2=wf2(i)
  t1=dble(z1)**2+aimag(z1)**2
  t2=dble(z2)**2+aimag(z2)**2
  z1=conjg(z1)*z2
  mag1(i)=mag1(i)+wo2*dble(z1)
  mag2(i)=mag2(i)+wo2*aimag(z1)
  mag3(i)=mag3(i)+wo*(t1-t2)
  rho(i)=rho(i)+wo*(t1+t2)
end do
return
end subroutine

pure subroutine rmk2(n,wo,wf1,wf2,rho,mag)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf1(n),wf2(n)
real(8), intent(inout) :: rho(n),mag(n)
! local variables
integer i
real(8) t1,t2
do i=1,n
  t1=dble(wf1(i))**2+aimag(wf1(i))**2
  t2=dble(wf2(i))**2+aimag(wf2(i))**2
  mag(i)=mag(i)+wo*(t1-t2)
  rho(i)=rho(i)+wo*(t1+t2)
end do
return
end subroutine

pure subroutine rmk3(n,wo,wf,rho)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: wo
complex(8), intent(in) :: wf(n)
real(8), intent(inout) :: rho(n)
rho(:)=rho(:)+wo*(dble(wf(:))**2+aimag(wf(:))**2)
return
end subroutine

end subroutine
!EOC

