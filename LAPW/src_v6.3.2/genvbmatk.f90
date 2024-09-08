
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvbmatk(vmt,vir,bmt,bir,ngp,igpig,wfmt,ld,wfgp,vbmat)
use modmain
use modomp
implicit none
! arguments
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
integer, intent(in) :: ld
complex(8), intent(in) :: wfgp(ld,nspinor,nstsv)
complex(8), intent(out) :: vbmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn
integer is,ias,nrc,nrci
integer npc,igp,nthd
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:),wfir(:,:),z(:)
! external functions
complex(8) zfcmtinp,zdotc
external zfcmtinp,zdotc
! zero the matrix elements
vbmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,ias,is,nrc) &
!$OMP PRIVATE(nrci,npc,ispn,ist) &
!$OMP NUM_THREADS(nthd)
allocate(wfmt1(npcmtmax,nspinor))
!$OMP DO
do jst=1,nstsv
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
      call vbmk1(npc,vmt(:,ias),bmt(:,ias,1),bmt(:,ias,2),bmt(:,ias,3),wfmt1, &
       wfmt1(:,2))
    else
! collinear case
      call vbmk2(npc,vmt(:,ias),bmt(:,ias,1),wfmt1,wfmt1(:,2))
    end if
    do ist=1,jst
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
call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir,z,ispn,jspn) &
!$OMP PRIVATE(igp,ist) &
!$OMP NUM_THREADS(nthd)
allocate(wfir(ngtot,nspinor),z(ngkmax))
!$OMP DO
do jst=1,nstsv
! Fourier transform wavefunction to real-space
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    wfir(:,ispn)=0.d0
    do igp=1,ngp(jspn)
      wfir(igfft(igpig(igp,jspn)),ispn)=wfgp(igp,ispn,jst)
    end do
    call zfftifc(3,ngridg,1,wfir(:,ispn))
  end do
! apply local potential and magnetic field to spinor wavefunction
  if (ncmag) then
! non-collinear case
    call vbmk1(ngtot,vir,bir,bir(:,2),bir(:,3),wfir,wfir(:,2))
  else
! collinear case
    call vbmk2(ngtot,vir,bir,wfir,wfir(:,2))
  end if
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! Fourier transform to G+p-space
    call zfftifc(3,ngridg,-1,wfir(:,ispn))
    do igp=1,ngp(jspn)
      z(igp)=wfir(igfft(igpig(igp,jspn)),ispn)
    end do
    do ist=1,jst
      vbmat(ist,jst)=vbmat(ist,jst)+zdotc(ngp(jspn),wfgp(:,ispn,ist),1,z,1)
    end do
  end do
end do
!$OMP END DO
deallocate(wfir,z)
!$OMP END PARALLEL
call freethd(nthd)
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vbmat(ist,jst)=conjg(vbmat(jst,ist))
  end do
end do
return

contains

pure subroutine vbmk1(n,v,b1,b2,b3,wf1,wf2)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: v(n),b1(n),b2(n),b3(n)
complex(8), intent(inout) :: wf1(n),wf2(n)
! local variables
integer i
real(8) t1
complex(8) z1,z2
do i=1,n
  z2=cmplx(b1(i),b2(i),8)
  t1=b3(i)
  z1=(v(i)+t1)*wf1(i)+conjg(z2)*wf2(i)
  wf2(i)=(v(i)-t1)*wf2(i)+z2*wf1(i)
  wf1(i)=z1
end do
return
end subroutine

pure subroutine vbmk2(n,v,b,wf1,wf2)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: v(n),b(n)
complex(8), intent(inout) :: wf1(n),wf2(n)
! local variables
integer i
do i=1,n
  wf1(i)=(v(i)+b(i))*wf1(i)
  wf2(i)=(v(i)-b(i))*wf2(i)
end do
return
end subroutine

end subroutine

