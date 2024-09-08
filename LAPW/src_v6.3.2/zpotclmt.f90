
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zpotclmt
! !INTERFACE:
subroutine zpotclmt(nr,nri,ld,rl,wpr,zrhomt,zvclmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr     : number of radial mesh points (in,integer)
!   nri    : number of points on inner part of muffin-tin (in,integer)
!   ld     : leading dimension (in,integer)
!   rl     : r^l on the radial mesh (in,real(ld,-lmaxo-1:lmaxo+2))
!   wpr    : weights for partial integration on radial mesh (in,real(4,nr))
!   zrhomt : muffin-tin charge density (in,complex(*))
!   zvclmt : muffin-tin Coulomb potential (out,complex(*))
! !DESCRIPTION:
!   Solves the Poisson equation for the charge density contained in an isolated
!   muffin-tin using the Green's function approach. In other words, the
!   spherical harmonic expansion of the Coulomb potential, $V_{lm}$, is obtained
!   from the density expansion, $\rho_{lm}$, by
!   $$ V_{lm}(r)=\frac{4\pi}{2l+1}\left(\frac{1}{r^{l+1}}\int_0^r\rho_{lm}(r')
!      {r'}^{l+2}dr'+r^l\int_r^R\frac{\rho_{lm}(r')}{{r'}^{l-1}}dr'\right) $$
!   where $R$ is the muffin-tin radius.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri,ld
real(8), intent(in) :: rl(ld,-lmaxo-1:lmaxo+2),wpr(4,nr)
complex(8), intent(in) :: zrhomt(*)
complex(8), intent(out) :: zvclmt(*)
! local variables
integer nro,iro,ir
integer l,l1,l2,l3
integer m,lm,npi,i
real(8) r1,r2,t0,t1,t2,t3,t4
! automatic arrays
real(8) f1(nr),f2(nr),f3(nr),f4(nr),f5(nr)
nro=nr-nri
iro=nri+1
npi=lmmaxi*nri
lm=0
do l=0,lmaxi
  l1=l+2
  l2=-l+1
  l3=-l-1
  t0=fourpi/dble(2*l+1)
  do m=-l,l
    lm=lm+1
    i=lm
    do ir=1,nri
      t1=dble(zrhomt(i)); t2=aimag(zrhomt(i))
      r1=rl(ir,l1); r2=rl(ir,l2)
      f1(ir)=t1*r1; f2(ir)=t2*r1
      f3(ir)=t1*r2; f4(ir)=t2*r2
      i=i+lmmaxi
    end do
    do ir=iro,nr
      t1=dble(zrhomt(i)); t2=aimag(zrhomt(i))
      r1=rl(ir,l1); r2=rl(ir,l2)
      f1(ir)=t1*r1; f2(ir)=t2*r1
      f3(ir)=t1*r2; f4(ir)=t2*r2
      i=i+lmmaxo
    end do
    call splintwp(nr,wpr,f1,f5)
    call splintwp(nr,wpr,f2,f1)
    call splintwp(nr,wpr,f3,f2)
    call splintwp(nr,wpr,f4,f3)
    t1=f2(nr); t2=f3(nr)
    i=lm
    do ir=1,nri
      r1=t0*rl(ir,l3); r2=t0*rl(ir,l)
      t3=r1*f5(ir)+r2*(t1-f2(ir))
      t4=r1*f1(ir)+r2*(t2-f3(ir))
      zvclmt(i)=cmplx(t3,t4,8)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      r1=t0*rl(ir,l3); r2=t0*rl(ir,l)
      t3=r1*f5(ir)+r2*(t1-f2(ir))
      t4=r1*f1(ir)+r2*(t2-f3(ir))
      zvclmt(i)=cmplx(t3,t4,8)
      i=i+lmmaxo
    end do
  end do
end do
do l=lmaxi+1,lmaxo
  l1=l+2
  l2=-l+1
  l3=-l-1
  t0=fourpi/dble(2*l+1)
  do m=-l,l
    lm=lm+1
    i=npi+lm
    do ir=iro,nr
      t1=dble(zrhomt(i)); t2=aimag(zrhomt(i))
      r1=rl(ir,l1); r2=rl(ir,l2)
      f1(ir)=t1*r1; f2(ir)=t2*r1
      f3(ir)=t1*r2; f4(ir)=t2*r2
      i=i+lmmaxo
    end do
    call splintwp(nro,wpr(1,iro),f1(iro),f5(iro))
    call splintwp(nro,wpr(1,iro),f2(iro),f1(iro))
    call splintwp(nro,wpr(1,iro),f3(iro),f2(iro))
    call splintwp(nro,wpr(1,iro),f4(iro),f3(iro))
    t1=f2(nr); t2=f3(nr)
    i=npi+lm
    do ir=iro,nr
      r1=t0*rl(ir,l3); r2=t0*rl(ir,l)
      t3=r1*f5(ir)+r2*(t1-f2(ir))
      t4=r1*f1(ir)+r2*(t2-f3(ir))
      zvclmt(i)=cmplx(t3,t4,8)
      i=i+lmmaxo
    end do
  end do
end do
return
end subroutine
!EOC

