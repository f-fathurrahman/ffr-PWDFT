
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzvbmatk(zvmt,zvir,zbmt,zbir,nst,ngp,igpig,wfmt,wfir,wfgp,vbmat)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(in) :: zvmt(npcmtmax,natmtot),zvir(ngtot)
complex(8), intent(in) :: zbmt(npcmtmax,natmtot,ndmag),zbir(ngtot,ndmag)
integer, intent(in) :: nst,ngp,igpig(ngkmax)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nst)
! note that wfir does not have a 1/sqrt(omega) prefactor
complex(8), intent(in) :: wfir(ngtot,nspinor,nst)
complex(8), intent(in) :: wfgp(ngkmax,nspinor,nst)
complex(8), intent(out) :: vbmat(nst,nst)
! local variables
integer ist,jst,ispn
integer is,ias,nrc,nrci
integer npc,nthd
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:),wfir1(:,:),z(:)
! external functions
complex(8) zfcmtinp,zdotc
external zfcmtinp,zdotc
! zero the matrix elements
vbmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
call holdthd(nst,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,ias,is,nrc) &
!$OMP PRIVATE(nrci,npc,ispn,ist) &
!$OMP NUM_THREADS(nthd)
allocate(wfmt1(npcmtmax,nspinor))
!$OMP DO
do jst=1,nst
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
    do ispn=1,nspinor
      call zcopy(npc,wfmt(:,ias,ispn,jst),1,wfmt1(:,ispn),1)
    end do
! apply local potential and magnetic field to spinor wavefunction
    if (ncmag) then
! non-collinear case
      call zvbmk1(npc,zvmt(:,ias),zbmt(:,ias,1),zbmt(:,ias,2),zbmt(:,ias,3), &
       wfmt1,wfmt1(:,2))
    else
! collinear case
      call zvbmk2(npc,zvmt(:,ias),zbmt(:,ias,1),wfmt1,wfmt1(:,2))
    end if
    do ist=1,nst
      do ispn=1,nspinor
! compute inner product (functions are in spherical coordinates)
        vbmat(ist,jst)=vbmat(ist,jst)+zfcmtinp(nrc,nrci,wrcmt(:,is), &
         wfmt(:,ias,ispn,ist),wfmt1(:,ispn))
      end do
    end do
  end do
end do
!$OMP END DO
deallocate(wfmt1)
!$OMP END PARALLEL
call freethd(nthd)
!---------------------------!
!     interstitial part     !
!---------------------------!
call holdthd(nst,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,z,ispn,ist) &
!$OMP NUM_THREADS(nthd)
allocate(wfir1(ngtot,nspinor),z(ngkmax))
!$OMP DO
do jst=1,nst
  do ispn=1,nspinor
    call zcopy(ngtot,wfir(:,ispn,jst),1,wfir1(:,ispn),1)
  end do
! apply local potential and magnetic field to spinor wavefunction
  if (ncmag) then
! non-collinear case
    call zvbmk1(ngtot,zvir,zbir,zbir(:,2),zbir(:,3),wfir1,wfir1(:,2))
  else
! collinear case
    call zvbmk2(ngtot,zvir,zbir,wfir1,wfir1(:,2))
  end if
  do ispn=1,nspinor
! Fourier transform to G+p-space
    call zfftifc(3,ngridg,-1,wfir1(:,ispn))
    z(1:ngp)=wfir1(igfft(igpig(1:ngp)),ispn)
    do ist=1,nst
      vbmat(ist,jst)=vbmat(ist,jst)+zdotc(ngp,wfgp(:,ispn,ist),1,z,1)
    end do
  end do
end do
!$OMP END DO
deallocate(wfir1,z)
!$OMP END PARALLEL
call freethd(nthd)
return

contains

pure subroutine zvbmk1(n,zv,zb1,zb2,zb3,wf1,wf2)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: zv(n),zb1(n),zb2(n),zb3(n)
complex(8), intent(inout) :: wf1(n),wf2(n)
! local variables
integer i
complex(8) v,b1,b2,z1
do i=1,n
  v=zv(i)
  b1=zb1(i)
  b2=cmplx(-aimag(zb2(i)),dble(zb2(i)),8)
  z1=(v+zb3(i))*wf1(i)+(b1-b2)*wf2(i)
  wf2(i)=(v-zb3(i))*wf2(i)+(b1+b2)*wf1(i)
  wf1(i)=z1
end do
return
end subroutine

pure subroutine zvbmk2(n,zv,zb,wf1,wf2)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: zv(n),zb(n)
complex(8), intent(inout) :: wf1(n),wf2(n)
! local variables
integer i
do i=1,n
  wf1(i)=(zv(i)+zb(i))*wf1(i)
  wf2(i)=(zv(i)-zb(i))*wf2(i)
end do
return
end subroutine

end subroutine

