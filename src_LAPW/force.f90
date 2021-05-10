SUBROUTINE force
  ! !USES:
  USE modmain
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
  IMPLICIT NONE 
  ! local variables
  INTEGER ik,idm,ispn
  INTEGER is,ias,nr,nri
  INTEGER np,i
  REAL(8) t1
  REAL(8) ts0,ts1
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: rfmt1(:,:),rfmt2(:),grfmt(:,:)
  ! external functions
  REAL(8) rfmtinp
  external rfmtinp

  CALL timesec(ts0)
  ALLOCATE(grfmt(npmtmax,3))
  !---------------------------------!
  !     Hellmann-Feynman forces     !
  !---------------------------------!
  DO ias=1,natmtot
    is=idxis(ias)
    nr=nrmt(is)
    nri=nrmti(is)
    CALL gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),vclmt(:,ias),npmtmax,grfmt)
    DO i=1,3
      forcehf(i,ias)=-spzn(is)*grfmt(1,i)*y00
    ENDDO 
  ENDDO 

  ! symmetrise Hellmann-Feynman forces
  CALL symveca(forcehf)
  
  !----------------------------------!
  !     IBS correction to forces     !
  !----------------------------------!
  ! set the IBS forces to zero
  forceibs(:,:)=0.d0
  ! compute k-point dependent contribution to the IBS forces
  DO ik=1,nkpt
    ! distribute among MPI processes
    CALL forcek(ik)
  ENDDO 

  ! integral of Kohn-Sham potential with gradient of density
  DO ias=1,natmtot
    is=idxis(ias)
    nr=nrmt(is)
    nri=nrmti(is)
    CALL gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rhomt(:,ias),npmtmax,grfmt)
    DO i=1,3
      t1=rfmtinp(nr,nri,wrmt(:,is),vsmt(:,ias),grfmt(:,i))
      forceibs(i,ias)=forceibs(i,ias)+t1
    ENDDO 
  ENDDO 

  ! integral of Kohn-Sham magnetic field with magnetisation gradient
  IF(spinpol) THEN 
    ALLOCATE(rfmt1(npmtmax,natmtot))
    DO idm=1,ndmag
      DO ias=1,natmtot
        is=idxis(ias)
        CALL rfsht(nrcmt(is),nrcmti(is),bsmt(:,ias,idm),rfmt1(:,ias))
      ENDDO 
      CALL rfmtctof(rfmt1)
      DO ias=1,natmtot
        is=idxis(ias)
        nr=nrmt(is)
        nri=nrmti(is)
        CALL gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),magmt(:,ias,idm), &
         npmtmax,grfmt)
        DO i=1,3
          t1=rfmtinp(nr,nri,wrmt(:,is),rfmt1(:,ias),grfmt(:,i))
          forceibs(i,ias)=forceibs(i,ias)+t1
        ENDDO 
      ENDDO 
    ENDDO 
    DEALLOCATE(rfmt1)
  ENDIF 

  ! integral of Kohn-Sham tau-DFT potential with kinetic energy density gradient
  IF(xcgrad.eq.4) THEN 
    ALLOCATE(rfmt1(npmtmax,natmtot),rfmt2(npmtmax))
    DO ias=1,natmtot
      is=idxis(ias)
      CALL rfsht(nrcmt(is),nrcmti(is),wsmt(:,ias),rfmt1(:,ias))
    ENDDO 
    CALL rfmtctof(rfmt1)
    DO ispn=1,nspinor
      DO ias=1,natmtot
        is=idxis(ias)
        nr=nrmt(is)
        nri=nrmti(is)
        np=npmt(is)
        rfmt2(1:np)=taumt(1:np,ias,ispn)-taucr(1:np,ias,ispn)
        CALL gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rfmt2,npmtmax,grfmt)
        DO i=1,3
          t1=rfmtinp(nr,nri,wrmt(:,is),rfmt1(:,ias),grfmt(:,i))
          forceibs(i,ias)=forceibs(i,ias)+t1
        ENDDO 
      ENDDO 
    ENDDO 
    DEALLOCATE(rfmt1,rfmt2)
  ENDIF 

  ! symmetrise IBS forces
  CALL symveca(forceibs)
  ! total force on each atom
  DO ias=1,natmtot
    forcetot(:,ias)=forcehf(:,ias)+forceibs(:,ias)
  ENDDO 
  ! symmetrise total forces
  CALL symveca(forcetot)

  ! compute the average force
  DO i=1,3
    forceav(i)=0.d0
    DO ias=1,natmtot
      forceav(i)=forceav(i)+forcetot(i,ias)
    ENDDO 
    forceav(i)=forceav(i)/dble(natmtot)
  ENDDO 

  ! remove the average force, if required, to prevent translation of atomic basis
  IF(tfav0) THEN 
    DO ias=1,natmtot
      forcetot(:,ias)=forcetot(:,ias)-forceav(:)
    ENDDO 
  ENDIF 

  ! zero force on atoms with negative mass
  DO ias=1,natmtot
    is=idxis(ias)
    IF(spmass(is).le.0.d0) forcetot(:,ias)=0.d0
  ENDDO 

  ! compute maximum force magnitude over all atoms
  forcemax=0.d0
  DO ias=1,natmtot
    t1=sqrt(forcetot(1,ias)**2+forcetot(2,ias)**2+forcetot(3,ias)**2)
    IF(t1.gt.forcemax) forcemax=t1
  ENDDO 
  DEALLOCATE(grfmt)

  CALL timesec(ts1)
  timefor=timefor+ts1-ts0
  RETURN 

END SUBROUTINE 