
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmatk(vmt,vir,ngp,igpig,wfmt,ld,wfgp,vmat)
use modmain

implicit none
! arguments
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
integer, intent(in) :: ld
complex(8), intent(in) :: wfgp(ld,nspinor,nstsv)
complex(8), intent(out) :: vmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn
integer is,ias,nrc,nrci
integer npc,igp
! allocatable arrays
complex(8), allocatable :: wfmt1(:),wfir(:),z(:)
! external functions
complex(8), external :: zfcmtinp,zdotc
! zero the matrix elements
vmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(wfmt1(npcmtmax))
do jst=1,nstsv
  do ispn=1,nspinor
    do ias=1,natmtot
      is=idxis(ias)
      nrc=nrcmt(is)
      nrci=nrcmti(is)
      npc=npcmt(is)
! apply potential to wavefunction
      wfmt1(1:npc)=vmt(1:npc,ias)*wfmt(1:npc,ias,ispn,jst)
      do ist=1,jst
! compute inner product (functions are in spherical coordinates)
        vmat(ist,jst)=vmat(ist,jst)+zfcmtinp(nrc,nrci,wrcmt(:,is), &
         wfmt(:,ias,ispn,ist),wfmt1)
      end do
    end do
  end do
end do
deallocate(wfmt1)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(wfir(ngtot),z(ngkmax))
do jst=1,nstsv
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! Fourier transform wavefunction to real-space
    wfir(:)=0.d0
    do igp=1,ngp(jspn)
      wfir(igfft(igpig(igp,jspn)))=wfgp(igp,ispn,jst)
    end do
    call zfftifc(3,ngridg,1,wfir)
! apply potential to wavefunction
    wfir(:)=vir(:)*wfir(:)
! Fourier transform to G+p-space
    call zfftifc(3,ngridg,-1,wfir)
    do igp=1,ngp(jspn)
      z(igp)=wfir(igfft(igpig(igp,jspn)))
    end do
    do ist=1,jst
! compute inner product
      vmat(ist,jst)=vmat(ist,jst)+zdotc(ngp(jspn),wfgp(:,ispn,ist),1,z,1)
    end do
  end do
end do
deallocate(wfir,z)
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vmat(ist,jst)=conjg(vmat(jst,ist))
  end do
end do
return
end subroutine

