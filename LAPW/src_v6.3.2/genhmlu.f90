
! Copyright (C) 2017 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genhmlu(ik0,h)
use modmain
use modulr
use modomp
implicit none
! arguments
integer, intent(in) :: ik0
complex(8), intent(out) :: h(nstulr,nstulr)
! local variables
integer ik,ist,jst,ispn,nthd
integer ikpa,jkpa,iq,ifq,igk
integer i1,i2,i3,j1,j2,j3,i,j
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:),wfgk(:,:,:)
complex(8), allocatable :: vmat(:,:)
! central k-point
ik=(ik0-1)*nkpa+1
! get the ground-state eigenvectors from file for central k-point
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
call getevecfv('.OUT',ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk(1,ik),vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
! calculate the wavefunctions for all states of the central k-point
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfgk(ngkmax,nspinor,nstsv))
call genwfsv(.false.,.true.,nstsv,idx,ngridg,igfft,ngk(:,ik),igkig(:,:,ik), &
 apwalm,evecfv,evecsv,wfmt,ngkmax,wfgk)
deallocate(apwalm,evecfv,evecsv)
! determine the interstitial wavefunctions in real-space (without 1/sqrt(omega))
allocate(wfir(ngtot,nspinor,nstsv))
call holdthd(nstsv,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ispn,igk) &
!$OMP NUM_THREADS(nthd)
do ist=1,nstsv
  do ispn=1,nspinor
    wfir(:,ispn,ist)=0.d0
    do igk=1,ngk(1,ik)
      wfir(igfft(igkig(igk,1,ik)),ispn,ist)=wfgk(igk,ispn,ist)
    end do
    call zfftifc(3,ngridg,1,wfir(:,ispn,ist))
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! generate the matrix elements for all Q-vectors
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vmat,iq,i,j,ikpa,jkpa) &
!$OMP PRIVATE(j1,j2,j3,ist,jst,i1,i2,i3) &
!$OMP NUM_THREADS(nthd)
allocate(vmat(nstsv,nstsv))
!$OMP DO
do ifq=1,nfqrz
  iq=iqrzf(ifq)
  if (spinpol) then
    call genzvbmatk(vsqmt(:,:,ifq),vsqir(:,ifq),bsqmt(:,:,:,ifq), &
     bsqir(:,:,ifq),nstsv,ngk(1,ik),igkig(:,1,ik),wfmt,wfir,wfgk,vmat)
  else
    call genzvmatk(vsqmt(:,:,ifq),vsqir(:,ifq),nstsv,ngk(1,ik),igkig(:,1,ik), &
     wfmt,wfir,wfgk,vmat)
  end if
  j=0
  do jkpa=1,nkpa
    j1=ivq(1,jkpa); j2=ivq(2,jkpa); j3=ivq(3,jkpa)
    do jst=1,nstsv
      j=j+1
      do ikpa=1,jkpa-1
        i=(ikpa-1)*nstsv+1
        i1=ivq(1,ikpa)-j1; i2=ivq(2,ikpa)-j2; i3=ivq(3,ikpa)-j3
        if (ivqiq(i1,i2,i3).eq.iq) then
          call zcopy(nstsv,vmat(:,jst),1,h(i,j),1)
        else if (ivqiq(-i1,-i2,-i3).eq.iq) then
          do ist=1,nstsv
            h(i,j)=conjg(vmat(jst,ist))
            i=i+1
          end do
        end if
      end do
      if (ifq.ne.1) cycle
      i=(jkpa-1)*nstsv+1
      call zcopy(jst,vmat(:,jst),1,h(i,j),1)
    end do
  end do
end do
!$OMP END DO
deallocate(vmat)
!$OMP END PARALLEL
call freethd(nthd)
! add the second-variational eigenvalues to the diagonal
i=0
do ikpa=1,nkpa
  ik=(ik0-1)*nkpa+ikpa
  do ist=1,nstsv
    i=i+1
    h(i,i)=h(i,i)+evalsv(ist,ik)
  end do
end do
deallocate(wfmt,wfir,wfgk)
return
end subroutine

