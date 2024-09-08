
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gencore
! !INTERFACE:
subroutine gencore
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Computes the core radial wavefunctions, eigenvalues and densities. The
!   radial Dirac equation is solved in the spherical part of the Kohn-Sham
!   potential to which the atomic potential has been appended for
!   $r>R_{\rm MT}$. In the case of spin-polarised calculations, and when
!   {\tt spincore} is set to {\tt .true.}, the Dirac equation is solved in the
!   spin-up and -down potentials created from the Kohn-Sham scalar potential and
!   magnetic field magnitude, with the occupancy divided equally between up and
!   down. The up and down densities determined in this way are added to both the
!   scalar density and the magnetisation in the routine {\tt rhocore}. Note
!   that this procedure is a simple, but inexact, approach to solving the radial
!   Dirac equation in a magnetic field.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Added polarised cores, November 2009 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ist,ispn,idm
integer is,ia,ja,ias,jas
integer nr,nri,nrs,ir,nthd
real(8) t1
! automatic arrays
logical done(natmmax)
real(8) vr(nrspmax),eval(nstspmax)
! allocatable arrays
real(8), allocatable :: br(:),fr(:,:)
if (spincore) allocate(br(nrmtmax),fr(nrmtmax,ndmag))
! loop over species and atoms
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  nrs=nrsp(is)
  done(:)=.false.
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
! Kohn-Sham magnetic field for spin-polarised core
    if (spincore) then
      do idm=1,ndmag
        call rfmtlm(1,nr,nri,bxcmt(:,ias,idm),fr(:,idm))
      end do
      if (ncmag) then
        do ir=1,nr
          br(ir)=sqrt(fr(ir,1)**2+fr(ir,2)**2+fr(ir,3)**2)*y00
        end do
      else
        do ir=1,nr
          br(ir)=abs(fr(ir,1))*y00
        end do
      end if
    end if
! loop over spin channels
    do ispn=1,nspncr
! use the spherical part of the crystal Kohn-Sham potential
      call rfmtlm(1,nr,nri,vsmt(:,ias),vr)
      vr(1:nr)=vr(1:nr)*y00
! spin-up and -down potentials for polarised core
      if (spincore) then
        if (ispn.eq.1) then
          vr(1:nr)=vr(1:nr)+br(1:nr)
        else
          vr(1:nr)=vr(1:nr)-br(1:nr)
        end if
      end if
! append the Kohn-Sham potential from the atomic calculation for r > R_MT
      t1=vr(nr)-vrsp(nr,is)
      do ir=nr+1,nrs
        vr(ir)=vrsp(ir,is)+t1
      end do
      rhocr(:,ias,ispn)=0.d0
      call holdthd(nstsp(is),nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(t1,ir) SHARED(is) &
!$OMP NUM_THREADS(nthd)
      do ist=1,nstsp(is)
        if (spcore(ist,is)) then
! solve the Dirac equation
          eval(ist)=evalcr(ist,ias)
          call rdirac(solsc,nsp(ist,is),lsp(ist,is),ksp(ist,is),nrs,rsp(:,is), &
           vr,eval(ist),rwfcr(:,1,ist,ias),rwfcr(:,2,ist,ias))
          if (spincore) then
! use the spin-averaged eigenvalue for the polarised core
            if (ispn.eq.1) then
              evalcr(ist,ias)=eval(ist)
            else
              evalcr(ist,ias)=0.5d0*(evalcr(ist,ias)+eval(ist))
            end if
            t1=0.5d0*occcr(ist,ias)
          else
            evalcr(ist,ias)=eval(ist)
            t1=occcr(ist,ias)
          end if
! add to the core density
!$OMP CRITICAL(gencore_)
          do ir=1,nr
            rhocr(ir,ias,ispn)=rhocr(ir,ias,ispn) &
             +t1*(rwfcr(ir,1,ist,ias)**2+rwfcr(ir,2,ist,ias)**2)
          end do
!$OMP END CRITICAL(gencore_)
        end if
      end do
!$OMP END PARALLEL DO
      call freethd(nthd)
      do ir=1,nr
        rhocr(ir,ias,ispn)=rhocr(ir,ias,ispn)*rlmt(ir,-2,is)*y00
      end do
! end loop over spin channels
    end do
    done(ia)=.true.
! copy to equivalent atoms
    do ja=1,natoms(is)
      if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
        jas=idxas(ja,is)
        do ist=1,nstsp(is)
          if (spcore(ist,is)) then
            evalcr(ist,jas)=evalcr(ist,ias)
            rwfcr(1:nrs,:,ist,jas)=rwfcr(1:nrs,:,ist,ias)
          end if
        end do
        rhocr(1:nr,jas,:)=rhocr(1:nr,ias,:)
        done(ja)=.true.
      end if
    end do
! end loop over species and atoms
  end do
end do
if (spincore) deallocate(br,fr)
return
end subroutine
!EOC

