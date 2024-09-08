
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genlofr
! !INTERFACE:
subroutine genlofr
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the local-orbital radial functions. This is done by integrating
!   the scalar relativistic Schr\"{o}dinger equation (or its energy deriatives)
!   at the current linearisation energies using the spherical part of the
!   Kohn-Sham potential. For each local-orbital, a linear combination of
!   {\tt lorbord} radial functions is constructed such that its radial
!   derivatives up to order ${\tt lorbord}-1$ are zero at the muffin-tin radius.
!   This function is normalised and the radial Hamiltonian applied to it. The
!   results are stored in the global array {\tt lofr}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!   Copied to equivalent atoms, February 2010 (A. Kozhevnikov and JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ja,ias,jas
integer nr,nri,ir,i
integer ilo,jlo,io,jo
integer nn,l,info
real(8) e,t1
! automatic arrays
logical done(natmmax)
integer ipiv(nplorb)
real(8) vr(nrmtmax),fr(nrmtmax)
real(8) p0(nrmtmax,lorbordmax),p1(nrmtmax)
real(8) q0(nrmtmax),q1(nrmtmax),ep0(nrmtmax,lorbordmax)
real(8) p0s(nrmtmax,nlomax),ep0s(nrmtmax,nlomax)
real(8) xa(nplorb),ya(nplorb)
real(8) a(nplorb,nplorb),b(nplorb)
! external functions
real(8) splint,polynm
external splint,polynm
do is=1,nspecies
  nr=nrmt(is)
  nri=nrmti(is)
  done(:)=.false.
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
! use spherical part of potential
    i=1
    do ir=1,nri
      vr(ir)=vsmt(i,ias)*y00
      i=i+lmmaxi
    end do
    do ir=nri+1,nr
      vr(ir)=vsmt(i,ias)*y00
      i=i+lmmaxo
    end do
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      do jo=1,lorbord(ilo,is)
! linearisation energy accounting for energy derivative
        e=lorbe(jo,ilo,ias)+dble(lorbdm(jo,ilo,is))*deapwlo
! integrate the radial Schrodinger equation
        call rschrodint(solsc,l,e,nr,rlmt(:,1,is),vr,nn,p0(:,jo),p1,q0,q1)
        ep0(1:nr,jo)=e*p0(1:nr,jo)
! normalise radial functions
        fr(1:nr)=p0(1:nr,jo)**2
        t1=splint(nr,rlmt(:,1,is),fr)
        t1=1.d0/sqrt(abs(t1))
        call dscal(nr,t1,p0(:,jo),1)
        call dscal(nr,t1,ep0(:,jo),1)
! set up the matrix of radial derivatives
        do i=1,nplorb
          ir=nr-nplorb+i
          xa(i)=rlmt(ir,1,is)
          ya(i)=p0(ir,jo)*rlmt(ir,-1,is)
        end do
        do io=1,lorbord(ilo,is)
          a(io,jo)=polynm(io-1,nplorb,xa,ya,rmt(is))
        end do
      end do
! set up the target vector
      b(:)=0.d0
      b(lorbord(ilo,is))=1.d0
      call dgesv(lorbord(ilo,is),1,a,nplorb,ipiv,b,nplorb,info)
      if (info.ne.0) goto 10
! generate linear superposition of radial functions
      p0s(:,ilo)=0.d0
      ep0s(:,ilo)=0.d0
      do io=1,lorbord(ilo,is)
        t1=b(io)
        call daxpy(nr,t1,p0(:,io),1,p0s(:,ilo),1)
        call daxpy(nr,t1,ep0(:,io),1,ep0s(:,ilo),1)
      end do
! normalise radial functions
      fr(1:nr)=p0s(1:nr,ilo)**2
      t1=splint(nr,rlmt(:,1,is),fr)
      t1=1.d0/sqrt(abs(t1))
      call dscal(nr,t1,p0s(:,ilo),1)
      call dscal(nr,t1,ep0s(:,ilo),1)
! subtract linear combination of previous local-orbitals with same l
      do jlo=1,ilo-1
        if (lorbl(jlo,is).eq.l) then
          fr(1:nr)=p0s(1:nr,ilo)*p0s(1:nr,jlo)
          t1=-splint(nr,rlmt(:,1,is),fr)
          call daxpy(nr,t1,p0s(:,jlo),1,p0s(:,ilo),1)
          call daxpy(nr,t1,ep0s(:,jlo),1,ep0s(:,ilo),1)
        end if
      end do
! normalise radial functions again
      fr(1:nr)=p0s(1:nr,ilo)**2
      t1=splint(nr,rlmt(:,1,is),fr)
      t1=abs(t1)
      if (t1.lt.1.d-25) goto 10
      t1=1.d0/sqrt(t1)
      call dscal(nr,t1,p0s(:,ilo),1)
      call dscal(nr,t1,ep0s(:,ilo),1)
! divide by r and store in global array
      do ir=1,nr
        t1=rlmt(ir,-1,is)
        lofr(ir,1,ilo,ias)=t1*p0s(ir,ilo)
        lofr(ir,2,ilo,ias)=t1*ep0s(ir,ilo)
      end do
    end do
    done(ia)=.true.
! copy to equivalent atoms
    do ja=1,natoms(is)
      if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
        jas=idxas(ja,is)
        do ilo=1,nlorb(is)
          call dcopy(nr,lofr(:,1,ilo,ias),1,lofr(:,1,ilo,jas),1)
          call dcopy(nr,lofr(:,2,ilo,ias),1,lofr(:,2,ilo,jas),1)
        end do
        done(ja)=.true.
      end if
    end do
! end loop over atoms and species
  end do
end do
return
10 continue
write(*,*)
write(*,'("Error(genlofr): degenerate local-orbital radial functions")')
write(*,'(" for species ",I4)') is
write(*,'(" atom ",I4)') ia
write(*,'(" and local-orbital ",I4)') ilo
write(*,*)
stop
end subroutine
!EOC

