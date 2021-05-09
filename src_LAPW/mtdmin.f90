! !DESCRIPTION:
!   Finds the atomic species pair for which the distance between the muffin-tin
!   surfaces is a minimum. This distance may be negative if the muffin-tins
!   overlap.
SUBROUTINE mtdmin(is,js,dmin)
  USE m_atoms, ONLY: atposc, nspecies, natoms
  USE m_lattice, ONLY: avec, epslat
  USE m_muffin_tins, ONLY: rmt
! !INPUT/OUTPUT PARAMETERS:
!   is, js : species numbers (out,integer)
!   dmin   : minimum distance between muffin-tin surfaces (out,real)
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(out) :: is,js
  REAL(8), intent(out) :: dmin
  ! local variables
  INTEGER :: i1,i2,i3,ks,ka,ls,la
  REAL(8) :: v1(3),v2(3),t1,t2,t3
  is=1
  js=1
  dmin=1.d6
  DO i1=-1,1
    DO i2=-1,1
      DO i3=-1,1
        v1(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3)
        DO ks=1,nspecies
          DO ka=1,natoms(ks)
            v2(:)=v1(:)+atposc(:,ka,ks)
            DO ls=1,nspecies
              t1=rmt(ks)+rmt(ls)
              DO la=1,natoms(ls)
                IF((i1.ne.0).or.(i2.ne.0).or.(i3.ne.0).or.(ks.ne.ls).or. &
                 (ka.ne.la)) THEN 
                  t2=sqrt((v2(1)-atposc(1,la,ls))**2 &
                         +(v2(2)-atposc(2,la,ls))**2 &
                         +(v2(3)-atposc(3,la,ls))**2)
                  t3=t2-t1
                  IF(t3.lt.dmin-epslat) THEN 
                    is=ks
                    js=ls
                    dmin=t3
                  ENDIF 
                ENDIF 
              ENDDO 
            ENDDO 
          ENDDO 
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 
