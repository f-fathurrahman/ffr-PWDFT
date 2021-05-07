
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepresk(ik,vclcv,vclvv,dvxmt,dvxir,dbxmt,dbxir)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: vclcv(ncrmax,natmtot,nstsv,nkpt)
complex(8), intent(in) :: vclvv(nstsv,nstsv,nkpt)
real(8), intent(inout) :: dvxmt(npcmtmax,natmtot),dvxir(ngtot)
real(8), intent(inout) :: dbxmt(npcmtmax,natmtot,ndmag),dbxir(ngtot,ndmag)
! local variables
integer ist,jst,idm
integer is,ia,ias,ic,m
integer nrc,nrci,npc
real(8) de
complex(8) z1,z2
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:),wfcr(:,:)
complex(8), allocatable :: zfmt1(:),zvfmt1(:,:)
complex(8), allocatable :: zfmt2(:,:),zfir2(:)
complex(8), allocatable :: zvfmt2(:,:,:),zvfir2(:,:)
! external functions
complex(8) rzfinp,rzfmtinp
external rzfinp,rzfmtinp
! get the eigenvalues/vectors from file for input k-point
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk(1,ik),vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
! calculate the wavefunctions for all states
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfir(ngtot,nspinor,nstsv))
call genwfsv(.false.,.false.,nstsv,idx,ngridg,igfft,ngk(1,ik),igkig(:,1,ik), &
 apwalm,evecfv,evecsv,wfmt,ngtot,wfir)
deallocate(apwalm,evecfv,evecsv)
!-----------------------------------------------------------!
!     core-conduction overlap density and magnetisation     !
!-----------------------------------------------------------!
allocate(wfcr(npcmtmax,2),zfmt1(npcmtmax))
if (spinpol) allocate(zvfmt1(npcmtmax,ndmag))
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    ic=0
    do ist=1,nstsp(is)
      if (spcore(ist,is)) then
        do m=-ksp(ist,is),ksp(ist,is)-1
          ic=ic+1
! pass in m-1/2 to wavefcr
          call wavefcr(.false.,lradstp,is,ia,ist,m,npcmtmax,wfcr)
          do jst=1,nstsv
            if (evalsv(jst,ik).gt.efermi) then
              if (spinpol) then
! compute the complex density and magnetisation
                call genzrm(npc,wfcr,wfcr(:,2),wfmt(:,ias,1,jst), &
                 wfmt(:,ias,2,jst),zfmt1,npcmtmax,zvfmt1)
              else
! compute the complex density
                zfmt1(1:npc)=conjg(wfcr(1:npc,1))*wfmt(1:npc,ias,1,jst)
              end if
              z1=conjg(vclcv(ic,ias,jst,ik))
              z2=rzfmtinp(nrc,nrci,wrcmt(:,is),vxmt(:,ias),zfmt1)
              z1=z1-conjg(z2)
              do idm=1,ndmag
                z2=rzfmtinp(nrc,nrci,wrcmt(:,is),bxmt(:,ias,idm),zvfmt1(:,idm))
                z1=z1-conjg(z2)
              end do
              de=evalcr(ist,ias)-evalsv(jst,ik)
              z1=z1*occmax*wkpt(ik)/(de+zi*swidth)
! residuals for exchange potential and field
!$OMP CRITICAL(oepresk_)
              call rzadd(npc,z1,zfmt1,dvxmt(:,ias))
              do idm=1,ndmag
                call rzadd(npc,z1,zvfmt1(:,idm),dbxmt(:,ias,idm))
              end do
!$OMP END CRITICAL(oepresk_)
! end loop over jst
            end if
          end do
        end do
! end loop over ist
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(wfcr,zfmt1)
if (spinpol) deallocate(zvfmt1)
!--------------------------------------------------------------!
!     valence-conduction overlap density and magnetisation     !
!--------------------------------------------------------------!
allocate(zfmt2(npcmtmax,natmtot),zfir2(ngtot))
if (spinpol) then
  allocate(zvfmt2(npcmtmax,natmtot,ndmag),zvfir2(ngtot,ndmag))
end if
do ist=1,nstsv
  if (evalsv(ist,ik).lt.efermi) then
    do jst=1,nstsv
      if (evalsv(jst,ik).gt.efermi) then
        if (spinpol) then
! compute the complex density and magnetisation
          call genzfrm(wfmt(:,:,1,ist),wfmt(:,:,2,ist),wfir(:,1,ist), &
           wfir(:,2,ist),wfmt(:,:,1,jst),wfmt(:,:,2,jst),wfir(:,1,jst), &
           wfir(:,2,jst),zfmt2,zfir2,zvfmt2,zvfir2)
        else
! compute the complex density
          call genzrho(.false.,.true.,ngtot,wfmt(:,:,:,ist),wfir(:,:,ist), &
           wfmt(:,:,:,jst),wfir(:,:,jst),zfmt2,zfir2)
        end if
        z1=conjg(vclvv(ist,jst,ik))
        z2=rzfinp(vxmt,vxir,zfmt2,zfir2)
        z1=z1-conjg(z2)
        do idm=1,ndmag
          z2=rzfinp(bxmt(:,:,idm),bxir(:,idm),zvfmt2(:,:,idm),zvfir2(:,idm))
          z1=z1-conjg(z2)
        end do
        de=evalsv(ist,ik)-evalsv(jst,ik)
        z1=z1*occmax*wkpt(ik)/(de+zi*swidth)
! add to residuals for exchange potential and field
!$OMP CRITICAL(oepresk_)
        call rzfadd(z1,zfmt2,zfir2,dvxmt,dvxir)
        do idm=1,ndmag
          call rzfadd(z1,zvfmt2(:,:,idm),zvfir2(:,idm),dbxmt(:,:,idm), &
           dbxir(:,idm))
        end do
!$OMP END CRITICAL(oepresk_)
! end loop over jst
      end if
    end do
! end loop over ist
  end if
end do
deallocate(wfmt,wfir,zfmt2,zfir2)
if (spinpol) deallocate(zvfmt2,zvfir2)
return
end subroutine

subroutine rzadd(n,za,zv,rv)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: za
real(8), intent(in) :: zv(2*n)
real(8), intent(out) :: rv(n)
! local variables
real(8) t1
t1=dble(za)
if (abs(t1).gt.1.d-12) call daxpy(n,t1,zv,2,rv,1)
t1=-aimag(za)
if (abs(t1).gt.1.d-12) call daxpy(n,t1,zv(2),2,rv,1)
return
end subroutine

subroutine rzfadd(za,zfmt,zfir,rfmt,rfir)
use modmain
implicit none
! arguments
complex(8), intent(in) :: za
complex(8), intent(in) :: zfmt(npcmtmax,natmtot),zfir(ngtot)
real(8), intent(inout) :: rfmt(npcmtmax,natmtot),rfir(ngtot)
! local variables
integer is,ias
do ias=1,natmtot
  is=idxis(ias)
  call rzadd(npcmt(is),za,zfmt(:,ias),rfmt(:,ias))
end do
call rzadd(ngtot,za,zfir,rfir)
return
end subroutine

