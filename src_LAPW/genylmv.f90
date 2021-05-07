! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   v    : input vector (in,real(3))
!   ylm  : array of spherical harmonics (out,complex((lmax+1)**2))
SUBROUTINE genylmv(lmax,v,ylm)
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: lmax
  REAL(8), intent(in) :: v(3)
  COMPLEX(8), intent(out) :: ylm(*)
  ! local variables
  INTEGER l,m,lm1,lm2,lm3,lm4
  REAL(8), parameter :: eps=1.d-14
  REAL(8) r,st,ct,sp,cp
  REAL(8) t1,t2,t3,t4
  COMPLEX(8) z1
  IF((lmax.lt.0).or.(lmax.gt.50)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(genylmv): lmax out of range : ",I8)') lmax
    WRITE(*,*)
    stop
  ENDIF 
  ylm(1)=0.28209479177387814347d0
  IF(lmax == 0) RETURN 
  
  r=sqrt(v(1)**2+v(2)**2+v(3)**2)
  IF(r > eps) THEN 
    t1 = v(3)/r
    IF(t1 >= 1.d0) THEN 
      st=0.d0
      ct=1.d0
    ELSEIF(t1 <= -1.d0) THEN 
      st=0.d0
      ct=-1.d0
    ELSE 
      st=sqrt(1.d0-t1**2)
      ct=t1
    ENDIF 
    IF((abs(v(1)) > eps) .or. (abs(v(2)) > eps)) THEN 
      t1=1.d0/sqrt(v(1)**2+v(2)**2)
      sp=t1*v(2)
      cp=t1*v(1)
    ELSE 
      sp=0.d0
      cp=1.d0
    ENDIF 
  ELSE 
    st=0.d0
    ct=1.d0
    sp=0.d0
    cp=1.d0
  ENDIF 
  z1=cmplx(cp,sp,8)
  ylm(3)=0.48860251190291992159d0*ct
  ylm(4)=-0.34549414947133547927d0*st*z1
  ylm(2)=-conjg(ylm(4))
  DO l=2,lmax
    lm1=(l+1)**2
    lm2=l**2
    lm3=(l-1)**2+1
    lm4=lm2+1
    ylm(lm1)=-st*sqrt(dble(2*l+1)/dble(2*l))*z1*ylm(lm2)
    IF(mod(l,2).eq.0) THEN 
      ylm(lm4)=conjg(ylm(lm1))
    ELSE 
      ylm(lm4)=-conjg(ylm(lm1))
    ENDIF 
    lm1=lm1-1
    ylm(lm1)=ct*sqrt(dble(2*l+1))*ylm(lm2)
    lm4=lm4+1
    IF(mod(l-1,2).eq.0) THEN 
      ylm(lm4)=conjg(ylm(lm1))
    ELSE 
      ylm(lm4)=-conjg(ylm(lm1))
    ENDIF 
    t1=ct*sqrt(dble((2*l-1)*(2*l+1)))
    t2=sqrt(dble((2*l+1))/dble(2*l-3))
    DO m=l-2,1,-1
      lm1=lm1-1
      lm2=lm2-1
      lm3=lm3-1
      lm4=lm4+1
      t3=1.d0/sqrt(dble((l-m)*(l+m)))
      t4=t2*sqrt(dble((l-m-1)*(l+m-1)))
      ylm(lm1)=t3*(t1*ylm(lm2)-t4*ylm(lm3))
      IF(mod(m,2).eq.0) THEN 
        ylm(lm4)=conjg(ylm(lm1))
      ELSE 
        ylm(lm4)=-conjg(ylm(lm1))
      ENDIF 
    ENDDO 
    lm1=lm1-1; lm2=lm2-1; lm3=lm3-1
    t3=1.d0/dble(l)
    t4=t2*dble(l-1)
    ylm(lm1)=t3*(t1*ylm(lm2)-t4*ylm(lm3))
  ENDDO 
  RETURN 
END SUBROUTINE 
