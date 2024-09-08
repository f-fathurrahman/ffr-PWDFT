
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rdirac
! !INTERFACE:
subroutine rdirac(sol,n,l,k,nr,r,vr,eval,g0,f0)
! !INPUT/OUTPUT PARAMETERS:
!   sol  : speed of light in atomic units (in,real)
!   n    : principal quantum number (in,integer)
!   l    : quantum number l (in,integer)
!   k    : quantum number k (l or l+1) (in,integer)
!   nr   : number of radial mesh points (in,integer)
!   r    : radial mesh (in,real(nr))
!   vr   : potential on radial mesh (in,real(nr))
!   eval : eigenvalue without rest-mass energy (inout,real)
!   g0   : major component of the radial wavefunction (out,real(nr))
!   f0   : minor component of the radial wavefunction (out,real(nr))
! !DESCRIPTION:
!   Finds the solution to the radial Dirac equation for a given potential $v(r)$
!   and quantum numbers $n$, $k$ and $l$. The method involves integrating the
!   equation using the predictor-corrector method and adjusting $E$ until the
!   number of nodes in the wavefunction equals $n-l-1$. The calling routine must
!   provide an initial estimate for the eigenvalue. Note that the arrays
!   {\tt g0} and {\tt f0} represent the radial functions multiplied by $r$.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: sol
integer, intent(in) :: n,l,k,nr
real(8), intent(in) :: r(nr),vr(nr)
real(8), intent(inout) :: eval
real(8), intent(out) :: g0(nr),f0(nr)
! local variables
integer, parameter :: maxit=2000
integer kpa,it,ir,irm
integer nn,nnd,nndp
! energy convergence tolerance
real(8), parameter :: eps=1.d-12
real(8) t1,de
! automatic arrays
real(8) g1(nr),f1(nr),fr(nr)
! external functions
real(8) splint
external splint
if (k.le.0) then
  write(*,*)
  write(*,'("Error(rdirac): k <= 0 : ",I8)') k
  write(*,*)
  stop
end if
if (k.gt.n) then
  write(*,*)
  write(*,'("Error(rdirac): incompatible n and k : ",2I8)') n,k
  write(*,*)
  stop
end if
if ((k.eq.n).and.(l.ne.k-1)) then
  write(*,*)
  write(*,'("Error(rdirac): incompatible n, k and l : ",3I8)') n,k,l
  write(*,*)
  stop
end if
if (k.eq.l) then
  kpa=k
else if (k.eq.l+1) then
  kpa=-k
else
  write(*,*)
  write(*,'("Error(rdirac): incompatible l and k : ",2I8)') l,k
  write(*,*)
  stop
end if
if (nr.lt.4) then
  write(*,*)
  write(*,'("Error(rdirac): nr < 4 : ",I8)') nr
  write(*,*)
  stop
end if
de=1.d0
nndp=0
do it=1,maxit
! integrate the Dirac equation
  call rdiracint(sol,kpa,eval,nr,r,vr,nn,g0,g1,f0,f1)
! check the number of nodes
  nnd=nn-(n-l-1)
  if (nnd.gt.0) then
    eval=eval-de
  else
    eval=eval+de
  end if
  if (it.gt.1) then
    if ((nnd.ne.0).or.(nndp.ne.0)) then
      if (nnd*nndp.le.0) then
        de=de*0.5d0
      else
        de=de*1.1d0
      end if
    end if
  end if
  nndp=nnd
  if (de.lt.eps*(abs(eval)+1.d0)) goto 20
end do
write(*,*)
write(*,'("Warning(rdirac): maximum iterations exceeded")')
20 continue
! find effective infinity and set wavefunction to zero after that point
! major component
irm=nr
do ir=2,nr
  if ((g0(ir-1)*g0(ir).lt.0.d0).or.(g1(ir-1)*g1(ir).lt.0.d0)) irm=ir
end do
g0(irm:nr)=0.d0
! minor component
irm=nr
do ir=2,nr
  if ((f0(ir-1)*f0(ir).lt.0.d0).or.(f1(ir-1)*f1(ir).lt.0.d0)) irm=ir
end do
f0(irm:nr)=0.d0
! normalise
do ir=1,nr
  fr(ir)=g0(ir)**2+f0(ir)**2
end do
t1=splint(nr,r,fr)
t1=sqrt(abs(t1))
if (t1.gt.0.d0) then
  t1=1.d0/t1
else
  write(*,*)
  write(*,'("Error(rdirac): zero wavefunction")')
  write(*,*)
  stop
end if
call dscal(nr,t1,g0,1)
call dscal(nr,t1,f0,1)
return
end subroutine
!EOC

