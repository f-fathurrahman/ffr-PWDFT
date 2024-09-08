
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmatk(tspndg,tlmdg,lmin,lmax,ias,ngp,apwalm,evecfv,evecsv,ld,dmat)
use modmain
implicit none
! arguments
logical, intent(in) :: tspndg,tlmdg
integer, intent(in) :: lmin,lmax
integer, intent(in) :: ias
integer, intent(in) :: ngp(nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
integer, intent(in) :: ld
complex(8), intent(out) :: dmat(ld,nspinor,ld,nspinor,nstsv)
! local variables
integer ist,ispn,jspn,is,ia
integer nrc,nrci,irco,irc
integer l,m1,m2,lm1,lm2
integer npc,npci,i1,i2,i,j
real(8) a,b,t1
complex(8) zq(2),z1
! automatic arrays
logical done(nstfv,nspnfv)
real(8) fr1(nrcmtmax),fr2(nrcmtmax)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:),wfmt2(:,:)
if (lmin.lt.0) then
  write(*,*)
  write(*,'("Error(gendmatk): lmin < 0 : ",I8)') lmin
  write(*,*)
  stop
end if
if (lmax.gt.lmaxo) then
  write(*,*)
  write(*,'("Error(gendmatk): lmax > lmaxo : ",2I8)') lmax,lmaxo
  write(*,*)
  stop
end if
! allocate local arrays
allocate(wfmt1(npcmtmax,nstfv,nspnfv),wfmt2(npcmtmax,nspinor))
! species and atom numbers
is=idxis(ias)
ia=idxia(ias)
nrc=nrcmt(is)
nrci=nrcmti(is)
irco=nrci+1
npc=npcmt(is)
npci=npcmti(is)
! de-phasing factor for spin-spirals
if (ssdph) then
  t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
  zq(1)=cmplx(cos(t1),sin(t1),8)
  zq(2)=conjg(zq(1))
end if
! zero the density matrix
dmat(:,:,:,:,:)=0.d0
done(:,:)=.false.
! begin loop over second-variational states
do j=1,nstsv
  if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
    wfmt2(1:npc,:)=0.d0
    i=0
    do ispn=1,nspinor
      jspn=jspnfv(ispn)
      do ist=1,nstfv
        i=i+1
        z1=evecsv(i,j)
        if (ssdph) z1=z1*zq(ispn)
        if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
          if (.not.done(ist,jspn)) then
            call wavefmt(lradstp,ias,ngp(jspn),apwalm(:,:,:,ias,jspn), &
             evecfv(:,ist,jspn),wfmt1(:,ist,jspn))
            done(ist,jspn)=.true.
          end if
! add to spinor wavefunction
          wfmt2(1:npc,ispn)=wfmt2(1:npc,ispn)+z1*wfmt1(1:npc,ist,jspn)
        end if
      end do
    end do
  else
! spin-unpolarised wavefunction
    call wavefmt(lradstp,ias,ngp,apwalm(:,:,:,ias,1),evecfv(:,j,1),wfmt2)
  end if
  do ispn=1,nspinor
    do jspn=1,nspinor
      if (tspndg.and.(ispn.ne.jspn)) cycle
      do l=lmin,lmax
        do m1=-l,l
          lm1=idxlm(l,m1)
          do m2=-l,l
            lm2=idxlm(l,m2)
            if (tlmdg.and.(lm1.ne.lm2)) cycle
            if (l.le.lmaxi) then
              i1=lm1; i2=lm2
              do irc=1,nrci
                z1=wfmt2(i1,ispn)*conjg(wfmt2(i2,jspn))
                fr1(irc)=dble(z1); fr2(irc)=aimag(z1)
                i1=i1+lmmaxi; i2=i2+lmmaxi
              end do
              do irc=irco,nrc
                z1=wfmt2(i1,ispn)*conjg(wfmt2(i2,jspn))
                fr1(irc)=dble(z1); fr2(irc)=aimag(z1)
                i1=i1+lmmaxo; i2=i2+lmmaxo
              end do
              a=dot_product(wrcmt(1:nrc,is),fr1(1:nrc))
              b=dot_product(wrcmt(1:nrc,is),fr2(1:nrc))
            else
              i1=npci+lm1; i2=npci+lm2
              do irc=irco,nrc
                z1=wfmt2(i1,ispn)*conjg(wfmt2(i2,jspn))
                fr1(irc)=dble(z1); fr2(irc)=aimag(z1)
                i1=i1+lmmaxo; i2=i2+lmmaxo
              end do
              a=dot_product(wrcmt(irco:nrc,is),fr1(irco:nrc))
              b=dot_product(wrcmt(irco:nrc,is),fr2(irco:nrc))
            end if
            dmat(lm1,ispn,lm2,jspn,j)=cmplx(a,b,8)
          end do
        end do
      end do
    end do
  end do
! end loop over second-variational states
end do
deallocate(wfmt1,wfmt2)
return
end subroutine

