! !INPUT/OUTPUT PARAMETERS:
!   nmatp  : order of overlap and Hamiltonian matrices (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   vpc    : p-vector in Cartesian coordinates (in,real(3))
!   h,o    : Hamiltonian and overlap matrices in packed or upper triangular
!            form (in,complex(*))
!   evalfv : first-variational eigenvalues (out,real(nstfv))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
! !DESCRIPTION:
!   This routine solves the first-variational eigenvalue equation for the
!   special case when inversion symmetry is present. In this case the
!   Hamiltonian and overlap matrices can be made real by using appropriate
!   linear combinations of the local-orbitals for atoms related by inversion
!   symmetry. These are derived from the effect of parity and complex
!   conjugation on the spherical harmonics: $P Y_{lm}=(-1)^l Y_{lm}$ and
!   $(Y_{lm})^*=(-1)^mY_{l-m}$.
SUBROUTINE eveqnfvr(nmatp,ngp,vpc,h,o,evalfv,evecfv)
  USE modmain
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nmatp,ngp
  REAL(8), intent(in) :: vpc(3)
  COMPLEX(8), intent(in) :: h(*),o(*)
  REAL(8), intent(out) :: evalfv(nstfv)
  COMPLEX(8), intent(out) :: evecfv(nmatmax,nstfv)
  ! local variables
  INTEGER :: is,ia,ja,jas
  INTEGER :: ilo,i,j,k,l,m
  INTEGER :: i1,i2,j1,j2
  INTEGER :: k1,k2,k3,k4
  INTEGER :: l1,l2,m1,m2
  INTEGER :: lwork,info
  REAL(8) :: v(3),vl,vu
  REAL(8) :: t1,t2,t3,t4
  REAL(8) :: ts0,ts1
  COMPLEX(8) :: h1,h2,o1,o2,z1

  ! ALLOCATABLE arrays
  LOGICAL, ALLOCATABLE :: tr(:),tp(:)
  INTEGER, ALLOCATABLE :: idx(:),s(:),map(:,:)
  INTEGER, ALLOCATABLE :: iwork(:),ifail(:)
  REAL(8), ALLOCATABLE :: rh(:),ro(:),w(:)
  REAL(8), ALLOCATABLE :: rv(:,:),work(:)
  COMPLEX(8), ALLOCATABLE :: zp(:)
  
  CALL timesec(ts0)
  ALLOCATE(tr(nlotot),tp(nlotot))
  ALLOCATE(idx(nlotot),s(nlotot))
  ALLOCATE(map(nlotot,nlotot))
  ALLOCATE(zp(nlotot))
  tp(:)=.false.
  i=0
  DO is=1,nspecies
    DO ia=1,natoms(is)
  ! symmetry equivalent atom, mapped with inversion
      ja=ieqatom(ia,is,2)
      jas=idxas(ja,is)
  ! residual phase factor
      v(:)=atposc(:,ia,is)+atposc(:,ja,is)
      t1=0.5d0*dot_product(vpc(:),v(:))
      z1=cmplx(cos(t1),sin(t1),8)
      DO ilo=1,nlorb(is)
        l=lorbl(ilo,is)
        DO m=-l,l
          i=i+1
  ! index to conjugate local-orbital in symmetry equivalent atom
          idx(i)=idxlo(idxlm(l,-m),ilo,jas)
          IF(ia.ne.ja) THEN 
  ! sign of parity and conjugation operators
            IF(mod(l+m,2).eq.0) THEN 
              s(i)=1
            ELSE 
              s(i)=-1
            ENDIF 
            IF(ia.lt.ja) THEN 
  ! if ia < ja use the real part of the sum of matrix elements
              tr(i)=.true.
            ELSEIF(ia.gt.ja) THEN 
  ! if ia > ja use the imaginary part of the difference of matrix elements
              s(i)=-s(i)
              tr(i)=.false.
            ENDIF 
          ELSE 
  ! if ia = ja THEN  use real function when l even and imaginary when l is odd
            IF(mod(m,2).eq.0) THEN 
              s(i)=1
            ELSE 
              s(i)=-1
            ENDIF 
  ! new function should be real if symmetric or imaginary if antisymmetric
            IF(mod(l,2).eq.0) THEN 
  ! l even
              IF(m.ge.0) THEN 
                tr(i)=.true.
              ELSE 
                s(i)=-s(i)
                tr(i)=.false.
              ENDIF 
            ELSE 
              ! l odd
              IF(m.ge.0) THEN 
                tr(i)=.false.
              ELSE 
                s(i)=-s(i)
                tr(i)=.true.
              ENDIF 
            ENDIF 
          ENDIF 
          ! phase factors if required
          IF(abs(t1) > 1.d-8) THEN 
            zp(i)=z1
            tp(i)=.true.
          ENDIF 
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO 
  ! map from local-orbital indices to position in matrix
  DO m=1,nlotot
    j=ngp+m
    DO l=1,m
      i=ngp+l
      map(l,m)=i+(j-1)*nmatp
      map(m,l)=map(l,m)
    ENDDO 
  ENDDO 
  !---------------------------------!
  !     real Hamiltonian matrix     !
  !---------------------------------!
  ALLOCATE(rh(nmatp**2))
  ! <APW|H|APW>
  DO j=1,ngp
    k=(j-1)*nmatp+1
    CALL dcopy(j,h(k),2,rh(k),1)
  ENDDO 
  ! <APW|H|lo>
  DO m1=1,nlotot
    j1=ngp+m1
    j2=ngp+idx(m1)
    DO i=1,ngp
      k1=i+(j1-1)*nmatp
      k2=i+(j2-1)*nmatp
      h1=h(k1); h2=h(k2)
      IF(tp(m1)) THEN 
        h1=h1*zp(m1); h2=h2*zp(m1)
      ENDIF 
      IF(tr(m1)) THEN 
        rh(k1)=dble(h1)+s(m1)*dble(h2)
      ELSE 
        rh(k1)=aimag(h1)+s(m1)*aimag(h2)
      ENDIF 
    ENDDO 
  ENDDO 
  ! <lo|H|lo>
  DO m1=1,nlotot
    m2=idx(m1)
    DO l1=1,m1
      l2=idx(l1)
      k1=map(l1,m1); k2=map(l1,m2); k3=map(l2,m1); k4=map(l2,m2)
      IF((tr(l1).and.tr(m1)).or.((.not.tr(l1)).and.(.not.tr(m1)))) THEN 
        rh(k1)=dble(h(k1))+s(m1)*dble(h(k2))+s(l1)*(dble(h(k3))+s(m1)*dble(h(k4)))
      ELSE 
        t2=aimag(h(k2))
        IF(l1.gt.m2) t2=-t2
        t3=aimag(h(k3))
        IF(l2.gt.m1) t3=-t3
        t4=aimag(h(k4))
        IF(l2.gt.m2) t4=-t4
        rh(k1)=aimag(h(k1))+s(m1)*t2+s(l1)*(t3+s(m1)*t4)
        IF(.not.tr(l1)) rh(k1)=-rh(k1)
      ENDIF 
    ENDDO 
  ENDDO 
  !-----------------------------!
  !     real overlap matrix     !
  !-----------------------------!
  ALLOCATE(ro(nmatp**2))
  ! <APW|O|APW>
  DO j=1,ngp
    k=(j-1)*nmatp+1
    CALL dcopy(j,o(k),2,ro(k),1)
  ENDDO 
  ! <APW|O|lo>
  DO m1=1,nlotot
    j1=ngp+m1
    j2=ngp+idx(m1)
    DO i=1,ngp
      k1=i+(j1-1)*nmatp
      k2=i+(j2-1)*nmatp
      o1=o(k1); o2=o(k2)
      IF(tp(m1)) THEN 
        o1=o1*zp(m1); o2=o2*zp(m1)
      ENDIF 
      IF(tr(m1)) THEN 
        ro(k1)=dble(o1)+s(m1)*dble(o2)
      ELSE 
        ro(k1)=aimag(o1)+s(m1)*aimag(o2)
      ENDIF 
    ENDDO 
  ENDDO 
  ! <lo|O|lo>
  DO m1=1,nlotot
    m2=idx(m1)
    DO l1=1,m1
      l2=idx(l1)
      k1=map(l1,m1); k2=map(l1,m2); k3=map(l2,m1); k4=map(l2,m2)
      IF((tr(l1).and.tr(m1)).or.((.not.tr(l1)).and.(.not.tr(m1)))) THEN 
        ro(k1)=dble(o(k1))+s(m1)*dble(o(k2))+s(l1)*(dble(o(k3))+s(m1)*dble(o(k4)))
      ELSE 
        t2=aimag(o(k2))
        IF(l1.gt.m2) t2=-t2
        t3=aimag(o(k3))
        IF(l2.gt.m1) t3=-t3
        t4=aimag(o(k4))
        IF(l2.gt.m2) t4=-t4
        ro(k1)=aimag(o(k1))+s(m1)*t2+s(l1)*(t3+s(m1)*t4)
        IF(.not.tr(l1)) ro(k1)=-ro(k1)
      ENDIF 
    ENDDO 
  ENDDO 
  ! solve the generalised eigenvalue problem for real symmetric matrices
  ALLOCATE(iwork(5*nmatp),ifail(nmatp))
  ALLOCATE(w(nmatp),rv(nmatp,nstfv))
  lwork=8*nmatp
  ALLOCATE(work(lwork))

  ! diagonalise the matrix
  CALL dsygvx(1,'V','I','U',nmatp,rh,nmatp,ro,nmatp,vl,vu,1,nstfv,evaltol,m,w, &
   rv,nmatp,work,lwork,iwork,ifail,info)

  IF(info.ne.0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(eveqnfvr): diagonalisation failed")')
    WRITE(*,'(" DSYGVX RETURN ed INFO = ",I8)') info
    IF(info.gt.nmatp) THEN 
      i=info-nmatp
      WRITE(*,'(" The leading minor of the overlap matrix of order ",I8)') i
      WRITE(*,'("  is not positive definite")')
      WRITE(*,'(" Order of overlap matrix : ",I8)') nmatp
    ENDIF 
    WRITE(*,*)
    STOP 
  ENDIF 
  evalfv(1:nstfv)=w(1:nstfv)
  ! reconstruct the complex eigenvectors
  DO j=1,nstfv
    evecfv(1:ngp,j)=rv(1:ngp,j)
    evecfv(ngp+1:nmatp,j)=0.d0
    DO l1=1,nlotot
      i1=ngp+l1
      i2=ngp+idx(l1)
      t1=rv(i1,j)
      IF(tr(l1)) THEN 
        evecfv(i1,j)=evecfv(i1,j)+t1
        evecfv(i2,j)=evecfv(i2,j)+s(l1)*t1
      ELSE 
        evecfv(i1,j)=evecfv(i1,j)-cmplx(0.d0,t1,8)
        evecfv(i2,j)=evecfv(i2,j)-cmplx(0.d0,s(l1)*t1,8)
      ENDIF 
    ENDDO 
    DO l1=1,nlotot
      IF(tp(l1)) THEN 
        i1=ngp+l1
        evecfv(i1,j)=evecfv(i1,j)*zp(l1)
      ENDIF 
    ENDDO 
  ENDDO 

  DEALLOCATE(iwork,ifail,w,rv,work)
  DEALLOCATE(tr,tp,idx,s,map,rh,ro,zp)

  CALL timesec(ts1)

  timefv=timefv+ts1-ts0
  RETURN 
END SUBROUTINE 
!EOC

