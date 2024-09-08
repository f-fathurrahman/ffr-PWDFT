
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: force
! !INTERFACE:
subroutine force
! !USES:
use modmain
use modtest
use modmpi
use modomp
! !DESCRIPTION:
!   Computes the various contributions to the atomic forces. In principle, the
!   force acting on a nucleus is simply the gradient at that site of the
!   classical electrostatic potential from the other nuclei and the electronic
!   density. This is a result of the Hellmann-Feynman theorem. However because
!   the basis set is dependent on the nuclear coordinates and is not complete,
!   the Hellman-Feynman force is inacurate and corrections to it are required.
!   The first is the core correction which arises because the core wavefunctions
!   were determined by neglecting the non-spherical parts of the Kohn-Sham
!   potential $v_s$. Explicitly this is given by
!   $$ {\bf F}_{\rm core}^{\alpha}=\int_{\rm MT_{\alpha}} v_s({\bf r})
!    \nabla\rho_{\rm core}^{\alpha}({\bf r})\,d{\bf r} $$
!   for atom $\alpha$. The second, which is the incomplete basis set (IBS)
!   correction, is due to the position dependence of the APW functions, and is
!   derived by considering the change in total energy if the eigenvector
!   coefficients were fixed and the APW functions themselves were changed. This
!   would result in changes to the first-variational Hamiltonian and overlap
!   matrices given by
!   \begin{align*}
!    \delta H_{\bf G,G'}^{\alpha}&=i({\bf G-G'})
!    \left(H^{\alpha}_{\bf G+k,G'+k}-\frac{1}{2}({\bf G+k})\cdot({\bf G'+k})
!    \tilde{\Theta}_{\alpha}({\bf G-G'})e^{-i({\bf G-G'})\cdot{\bf r}_{\alpha}}
!    \right)\\
!    \delta O_{\bf G,G'}^{\alpha}&=i({\bf G-G'})\left(O^{\alpha}_{\bf G+k,G'+k}
!    -\tilde{\Theta}_{\alpha}({\bf G-G'})e^{-i({\bf G-G'})\cdot{\bf r}_{\alpha}}
!    \right)
!   \end{align*}
!   where both ${\bf G}$ and ${\bf G'}$ run over the APW indices;
!   $\tilde{\Theta}_{\alpha}$ is the form factor of the smooth step function for
!   muffin-tin $\alpha$; and $H^{\alpha}$ and $O^{\alpha}$ are the muffin-tin
!   Hamiltonian and overlap matrices, respectively. The APW-local-orbital part
!   is given by
!   \begin{align*}
!    \delta H_{\bf G,G'}^{\alpha}&=i({\bf G+k})H^{\alpha}_{\bf G+k,G'+k}\\
!    \delta O_{\bf G,G'}^{\alpha}&=i({\bf G+k})O^{\alpha}_{\bf G+k,G'+k}
!   \end{align*}
!   where ${\bf G}$ runs over the APW indices and ${\bf G'}$ runs over the
!   local-orbital indices. There is no contribution from the
!   local-orbital-local-orbital part of the matrices. We can now write the IBS
!   correction in terms of the basis of first-variational states as
!   \begin{align*}
!    {\bf F}_{ij}^{\alpha{\bf k}}=\sum_{\bf G,G'}
!    b^{i{\bf k}*}_{\bf G}b^{j{\bf k}}_{\bf G'}\left(
!    \delta H_{\bf G,G'}^{\alpha}-\epsilon_j\delta O_{\bf G,G'}^{\alpha}\right),
!   \end{align*}
!   where $b^{i{\bf k}}$ is the first-variational eigenvector.
!   Finally, the ${\bf F}_{ij}^{\alpha{\bf k}}$ matrix elements can be
!   multiplied by the second-variational coefficients, and contracted over all
!   indices to obtain the IBS force:
!   \begin{align*}
!    {\bf F}_{\rm IBS}^{\alpha}=\sum_{\bf k}w_{\bf k}\sum_{l\sigma}n_{l{\bf k}}
!    \sum_{ij}c_{\sigma i}^{l{\bf k}*}c_{\sigma j}^{l{\bf k}}
!    {\bf F}_{ij}^{\alpha{\bf k}}
!    +\int_{\rm MT_{\alpha}}v_s({\bf r})\nabla\left[\rho({\bf r})
!    -\rho^{\alpha}_{\rm core}({\bf r})\right]\,d{\bf r},
!   \end{align*}
!   where $c^{l{\bf k}}$ are the second-variational coefficients, $w_{\bf k}$
!   are the $k$-point weights, $n_{l{\bf k}}$ are the occupancies. See routines
!   {\tt hmlaa}, {\tt olpaa}, {\tt hmlalo}, {\tt olpalo}, {\tt energy},
!   {\tt eveqn} and {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!   Fixed problem with second-variational forces, May 2008 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,idm,ispn
integer is,ias,nr,nri
integer np,i,nthd
real(8) t1
real(8) ts0,ts1
! allocatable arrays
real(8), allocatable :: rfmt1(:,:),rfmt2(:),grfmt(:,:)
! external functions
real(8) rfmtinp
external rfmtinp
call timesec(ts0)
allocate(grfmt(npmtmax,3))
!---------------------------------!
!     Hellmann-Feynman forces     !
!---------------------------------!
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),vclmt(:,ias),npmtmax,grfmt)
  do i=1,3
    forcehf(i,ias)=-spzn(is)*grfmt(1,i)*y00
  end do
end do
! symmetrise Hellmann-Feynman forces
call symveca(forcehf)
!----------------------------------!
!     IBS correction to forces     !
!----------------------------------!
! set the IBS forces to zero
forceibs(:,:)=0.d0
! compute k-point dependent contribution to the IBS forces
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  call forcek(ik)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! add IBS forces from each process and redistribute
if (np_mpi.gt.1) then
  call mpi_allreduce(mpi_in_place,forceibs,3*natmtot,mpi_double_precision, &
   mpi_sum,mpicom,ierror)
end if
! integral of Kohn-Sham potential with gradient of density
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rhomt(:,ias),npmtmax,grfmt)
  do i=1,3
    t1=rfmtinp(nr,nri,wrmt(:,is),vsmt(:,ias),grfmt(:,i))
    forceibs(i,ias)=forceibs(i,ias)+t1
  end do
end do
! integral of Kohn-Sham magnetic field with magnetisation gradient
if (spinpol) then
  allocate(rfmt1(npmtmax,natmtot))
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      call rfsht(nrcmt(is),nrcmti(is),bsmt(:,ias,idm),rfmt1(:,ias))
    end do
    call rfmtctof(rfmt1)
    do ias=1,natmtot
      is=idxis(ias)
      nr=nrmt(is)
      nri=nrmti(is)
      call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),magmt(:,ias,idm), &
       npmtmax,grfmt)
      do i=1,3
        t1=rfmtinp(nr,nri,wrmt(:,is),rfmt1(:,ias),grfmt(:,i))
        forceibs(i,ias)=forceibs(i,ias)+t1
      end do
    end do
  end do
  deallocate(rfmt1)
end if
! integral of Kohn-Sham tau-DFT potential with kinetic energy density gradient
if (xcgrad.eq.4) then
  allocate(rfmt1(npmtmax,natmtot),rfmt2(npmtmax))
  do ias=1,natmtot
    is=idxis(ias)
    call rfsht(nrcmt(is),nrcmti(is),wsmt(:,ias),rfmt1(:,ias))
  end do
  call rfmtctof(rfmt1)
  do ispn=1,nspinor
    do ias=1,natmtot
      is=idxis(ias)
      nr=nrmt(is)
      nri=nrmti(is)
      np=npmt(is)
      rfmt2(1:np)=taumt(1:np,ias,ispn)-taucr(1:np,ias,ispn)
      call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rfmt2,npmtmax,grfmt)
      do i=1,3
        t1=rfmtinp(nr,nri,wrmt(:,is),rfmt1(:,ias),grfmt(:,i))
        forceibs(i,ias)=forceibs(i,ias)+t1
      end do
    end do
  end do
  deallocate(rfmt1,rfmt2)
end if
! symmetrise IBS forces
call symveca(forceibs)
! total force on each atom
do ias=1,natmtot
  forcetot(:,ias)=forcehf(:,ias)+forceibs(:,ias)
end do
! symmetrise total forces
call symveca(forcetot)
! compute the average force
do i=1,3
  forceav(i)=0.d0
  do ias=1,natmtot
    forceav(i)=forceav(i)+forcetot(i,ias)
  end do
  forceav(i)=forceav(i)/dble(natmtot)
end do
! remove the average force, if required, to prevent translation of atomic basis
if (tfav0) then
  do ias=1,natmtot
    forcetot(:,ias)=forcetot(:,ias)-forceav(:)
  end do
end if
! zero force on atoms with negative mass
do ias=1,natmtot
  is=idxis(ias)
  if (spmass(is).le.0.d0) forcetot(:,ias)=0.d0
end do
! compute maximum force magnitude over all atoms
forcemax=0.d0
do ias=1,natmtot
  t1=sqrt(forcetot(1,ias)**2+forcetot(2,ias)**2+forcetot(3,ias)**2)
  if (t1.gt.forcemax) forcemax=t1
end do
deallocate(grfmt)
call timesec(ts1)
timefor=timefor+ts1-ts0
! write total forces to test file
call writetest(750,'total forces',nv=3*natmtot,tol=1.d-3,rva=forcetot)
return
end subroutine
!EOC

