SUBROUTINE plotpt3d(vpl)
  USE m_plotting, ONLY: vclp3d, np3d
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(out) :: vpl(3,np3d(1)*np3d(2)*np3d(3))
  ! local variables
  INTEGER ip,i1,i2,i3
  REAL(8) v1(3),v2(3),v3(3)
  REAL(8) t1,t2,t3
  ! generate 3D grid from corner vectors
  v1(:)=vclp3d(:,1)-vclp3d(:,0)
  v2(:)=vclp3d(:,2)-vclp3d(:,0)
  v3(:)=vclp3d(:,3)-vclp3d(:,0)
  ip=0
  DO i3=0,np3d(3)-1
    t3=dble(i3)/dble(np3d(3))
    DO i2=0,np3d(2)-1
      t2=dble(i2)/dble(np3d(2))
      DO i1=0,np3d(1)-1
        t1=dble(i1)/dble(np3d(1))
        ip=ip+1
        vpl(:,ip)=t1*v1(:)+t2*v2(:)+t3*v3(:)+vclp3d(:,0)
      ENDDO 
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 