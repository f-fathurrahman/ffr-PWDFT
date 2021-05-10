SUBROUTINE forcek(ik)
  ! !USES:
  use modmain
  ! !INPUT/OUTPUT PARAMETERS:
  !   ik : reduced k-point number (in,integer)
  ! !DESCRIPTION:
  !   Computes the {\bf k}-dependent contribution to the incomplete basis set
  !   (IBS) force. See the calling routine {\tt force} for a full description.
IMPLICIT NONE 
! arguments
INTEGER, intent(in) :: ik
! local variables
INTEGER ispn0,ispn1,ispn,jspn
INTEGER n,nm,nm2,is,ias,ist,jst
INTEGER iv(3),jv(3),ig,i,j,k,l
REAL(8) vj(3),sum,t1
COMPLEX(8) z1,z2
! ALLOCATABLE arrays
INTEGER, ALLOCATABLE :: ijg(:)
REAL(8), ALLOCATABLE :: dp(:),evalfv(:,:)
COMPLEX(8), ALLOCATABLE :: apwalm(:,:,:,:)
COMPLEX(8), ALLOCATABLE :: evecfv(:,:,:),evecsv(:,:)
COMPLEX(8), ALLOCATABLE :: h(:),o(:),dlh(:),dlo(:)
COMPLEX(8), ALLOCATABLE :: vh(:),vo(:),ffv(:,:),y(:)
! external functions
COMPLEX(8) zdotc
external zdotc
nm2=nmatmax**2
! allocate local arrays
ALLOCATE(ijg(nm2),dp(nm2))
ALLOCATE(evalfv(nstfv,nspnfv))
ALLOCATE(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
ALLOCATE(evecfv(nmatmax,nstfv,nspnfv))
ALLOCATE(h(nm2),o(nm2),dlh(nm2),dlo(nm2))
ALLOCATE(vh(nmatmax),vo(nmatmax))
ALLOCATE(ffv(nstfv,nstfv),y(nstfv))

! get the eigenvalues/vectors from file
CALL getevalfv(filext,ik,vkl(:,ik),evalfv)
CALL getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
IF(tevecsv) THEN 
  ALLOCATE(evecsv(nstsv,nstsv))
  CALL getevecsv(filext,ik,vkl(:,ik),evecsv)
ENDIF 

! loop over first-variational spin components
DO jspn=1,nspnfv
  IF(spinsprl) THEN 
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  ENDIF 
  n=ngk(jspn,ik)
  nm=nmat(jspn,ik)
  DO j=1,n
    k=(j-1)*nm
    jv(:)=ivg(:,igkig(j,jspn,ik))
    vj(:)=0.5d0*vgkc(:,j,jspn,ik)
    DO i=1,j
      k=k+1
      iv(:)=ivg(:,igkig(i,jspn,ik))-jv(:)
      ijg(k)=ivgig(iv(1),iv(2),iv(3))
      dp(k)=dot_product(vgkc(:,i,jspn,ik),vj(:))
    ENDDO 
  ENDDO 
! find the matching coefficients
  CALL match(n,vgkc(:,:,jspn,ik),gkc(:,jspn,ik),sfacgk(:,:,jspn,ik),apwalm)
! loop over species and atoms
  DO ias=1,natmtot
    is=idxis(ias)
! Hamiltonian and overlap matrices
    h(:)=0.d0
    CALL hmlaa(.false.,ias,n,apwalm(:,:,:,ias),nm,h)
    CALL hmlalo(ias,n,apwalm(:,:,:,ias),nm,h)
    o(:)=0.d0
    CALL olpaa(.false.,ias,n,apwalm(:,:,:,ias),nm,o)
    CALL olpalo(ias,n,apwalm(:,:,:,ias),nm,o)
! loop over Cartesian directions
    DO l=1,3
! APW-APW contribution
      DO j=1,n
        k=(j-1)*nm
        DO i=1,j
          k=k+1
          ig=ijg(k)
          t1=vgc(l,ig)
          z1=-ffacg(ig,is)*conjg(sfacg(ig,ias))
          z2=t1*(dp(k)*z1+h(k))
          dlh(k)=cmplx(-aimag(z2),dble(z2),8)
          z2=t1*(z1+o(k))
          dlo(k)=cmplx(-aimag(z2),dble(z2),8)
        ENDDO 
      ENDDO 
      DO j=n+1,nm
        k=(j-1)*nm
! APW-local-orbital contribution
        DO i=1,n
          k=k+1
          t1=vgkc(l,i,jspn,ik)
          z1=t1*h(k)
          dlh(k)=cmplx(-aimag(z1),dble(z1),8)
          z1=t1*o(k)
          dlo(k)=cmplx(-aimag(z1),dble(z1),8)
        ENDDO 
! zero the local-orbital-local-orbital contribution
        DO i=n+1,j
          k=k+1
          dlh(k)=0.d0
          dlo(k)=0.d0
        ENDDO 
      ENDDO 
! compute the force matrix elements in the first-variational basis
      DO jst=1,nstfv
        CALL zhemv('U',nm,zone,dlh,nm,evecfv(:,jst,jspn),1,zzero,vh,1)
        CALL zhemv('U',nm,zone,dlo,nm,evecfv(:,jst,jspn),1,zzero,vo,1)
        t1=evalfv(jst,jspn)
        DO ist=1,nstfv
          z1=zdotc(nm,evecfv(:,ist,jspn),1,vh,1)
          z2=zdotc(nm,evecfv(:,ist,jspn),1,vo,1)
          ffv(ist,jst)=z1-t1*z2
        ENDDO 
      ENDDO 
! compute the force using the second-variational coefficients if required
      sum=0.d0
      IF(tevecsv) THEN 
! spin-polarised case
        DO j=1,nstsv
          DO ispn=ispn0,ispn1
            i=(ispn-1)*nstfv+1
            CALL zgemv('N',nstfv,nstfv,zone,ffv,nstfv,evecsv(i,j),1,zzero,y,1)
            z1=zdotc(nstfv,evecsv(i,j),1,y,1)
            sum=sum+occsv(j,ik)*dble(z1)
          ENDDO 
        ENDDO 
      ELSE 
        ! spin-unpolarised case
        DO j=1,nstsv
          sum=sum+occsv(j,ik)*dble(ffv(j,j))
        ENDDO 
      ENDIF 

      forceibs(l,ias)=forceibs(l,ias)+wkpt(ik)*sum
! end loop over Cartesian components
    ENDDO 
! end loop over atoms and species
  ENDDO 
! end loop over first-variational spins
ENDDO 
DEALLOCATE(ijg,dp,evalfv,apwalm,evecfv)
IF(tevecsv) DEALLOCATE(evecsv)
DEALLOCATE(h,o,dlh,dlo,vh,vo,ffv,y)
RETURN 
END SUBROUTINE 
