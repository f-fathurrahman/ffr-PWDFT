SUBROUTINE findsymlat()
  USE m_symmetry, ONLY: isymlat, symlat, symlatd, nsymlat, symlatc
  USE m_spin, ONLY: vqlss, spinsprl
  USE m_lattice, ONLY: epslat, ainv, avec
  USE m_electric_vector_pot, ONLY: efieldl, afieldl, tefield, tafield
  ! !DESCRIPTION:
  !   Finds the point group symmetries which leave the Bravais lattice invariant.
  !   Let $A$ be the matrix consisting of the lattice vectors in columns, THEN 
  !   $$ g=A^{\rm T}A $$
  !   is the metric tensor. Any $3\times 3$ matrix $S$ with elements $-1$, 0 or 1
  !   is a point group symmetry of the lattice if $\det(S)$ is $-1$ or 1, and
  !   $$ S^{\rm T}gS=g. $$
  !   The first matrix in the set RETURN ed is the identity.
  IMPLICIT NONE 
  ! local variables
  INTEGER :: md,sym(3,3),i,j
  INTEGER :: i11,i12,i13,i21,i22,i23,i31,i32,i33
  REAL(8) :: s(3,3),g(3,3),sgs(3,3)
  REAL(8) :: c(3,3),v(3),t1
  ! external functions
  INTEGER i3mdet
  external i3mdet
  ! determine metric tensor
  CALL r3mtm(avec,avec,g)
  ! loop over all possible symmetry matrices
  nsymlat=0
  DO i11=-1,1; DO i12=-1,1; DO i13=-1,1
  DO i21=-1,1; DO i22=-1,1; DO i23=-1,1
  DO i31=-1,1; DO i32=-1,1; DO i33=-1,1
    sym(1,1)=i11; sym(1,2)=i12; sym(1,3)=i13
    sym(2,1)=i21; sym(2,2)=i22; sym(2,3)=i23
    sym(3,1)=i31; sym(3,2)=i32; sym(3,3)=i33
  ! determinant of matrix
    md=i3mdet(sym)
  ! matrix should be orthogonal
    IF(abs(md).ne.1) goto 10
  ! check invariance of metric tensor
    s(:,:)=dble(sym(:,:))
    CALL r3mtm(s,g,c)
    CALL r3mm(c,s,sgs)
    DO j=1,3
      DO i=1,3
        IF(abs(sgs(i,j)-g(i,j)).gt.epslat) goto 10
      ENDDO 
    ENDDO 
  ! check invariance of spin-spiral q-vector if required
    IF(spinsprl) THEN 
      CALL r3mtv(s,vqlss,v)
      t1=abs(vqlss(1)-v(1))+abs(vqlss(2)-v(2))+abs(vqlss(3)-v(3))
      IF(t1.gt.epslat) goto 10
    ENDIF 
  ! check invariance of electric field if required
    IF(tefield) THEN 
      CALL r3mv(s,efieldl,v)
      t1=abs(efieldl(1)-v(1))+abs(efieldl(2)-v(2))+abs(efieldl(3)-v(3))
      IF(t1.gt.epslat) goto 10
    ENDIF 
    
    ! check invariance of A-field if required
    IF(tafield) THEN 
      CALL r3mv(s,afieldl,v)
      t1=abs(afieldl(1)-v(1))+abs(afieldl(2)-v(2))+abs(afieldl(3)-v(3))
      IF(t1.gt.epslat) goto 10
    ENDIF 
    
    nsymlat=nsymlat+1
    IF(nsymlat > 48) THEN 
      WRITE(*,*)
      WRITE(*,'("Error(findsymlat): more than 48 symmetries found")')
      WRITE(*,'(" (lattice vectors may be linearly dependent)")')
      WRITE(*,*)
      stop
    ENDIF 
  ! store the symmetry and determinant in global arrays
    symlat(:,:,nsymlat)=sym(:,:)
    symlatd(nsymlat)=md
  10 continue
  ENDDO ; ENDDO ; ENDDO 
  ENDDO ; ENDDO ; ENDDO 
  ENDDO ; ENDDO ; ENDDO 
  IF(nsymlat.eq.0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(findsymlat): no symmetries found")')
    WRITE(*,*)
    stop
  ENDIF 
  ! make the first symmetry the identity
  DO i=1,nsymlat
    IF((symlat(1,1,i).eq.1).and.(symlat(1,2,i).eq.0).and.(symlat(1,3,i).eq.0) &
   .and.(symlat(2,1,i).eq.0).and.(symlat(2,2,i).eq.1).and.(symlat(2,3,i).eq.0) &
   .and.(symlat(3,1,i).eq.0).and.(symlat(3,2,i).eq.0).and.(symlat(3,3,i).eq.1)) &
    THEN 
      sym(:,:)=symlat(:,:,1)
      symlat(:,:,1)=symlat(:,:,i)
      symlat(:,:,i)=sym(:,:)
      md=symlatd(1)
      symlatd(1)=symlatd(i)
      symlatd(i)=md
      exit
    ENDIF 
  ENDDO 
  ! index to the inverse of each operation
  DO i=1,nsymlat
    CALL i3minv(symlat(:,:,i),sym)
    DO j=1,nsymlat
      IF((symlat(1,1,j).eq.sym(1,1)).and.(symlat(1,2,j).eq.sym(1,2)).and. &
          (symlat(1,3,j).eq.sym(1,3)).and.(symlat(2,1,j).eq.sym(2,1)).and. &
          (symlat(2,2,j).eq.sym(2,2)).and.(symlat(2,3,j).eq.sym(2,3)).and. &
          (symlat(3,1,j).eq.sym(3,1)).and.(symlat(3,2,j).eq.sym(3,2)).and. &
          (symlat(3,3,j).eq.sym(3,3))) THEN 
        isymlat(i)=j
        goto 30
      ENDIF 
    ENDDO 
    WRITE(*,*)
    WRITE(*,'("Error(findsymlat): inverse operation not found")')
    WRITE(*,'(" for lattice symmetry ",I2)') i
    WRITE(*,*)
    stop
  30 continue
  ENDDO 
  
  ! determine the lattice symmetries in Cartesian coordinates
  DO i=1,nsymlat
    s(:,:)=dble(symlat(:,:,i))
    CALL r3mm(s,ainv,c)
    CALL r3mm(avec,c,symlatc(:,:,i))
  ENDDO 
  RETURN 
END SUBROUTINE 
