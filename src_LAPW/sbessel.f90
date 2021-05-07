SUBROUTINE sbessel(lmax,x,jl)
  ! !INPUT/OUTPUT PARAMETERS:
  !   lmax : maximum order of Bessel function (in,integer)
  !   x    : real argument (in,real)
  !   jl   : array of RETURN ed values (out,real(0:lmax))
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: lmax
  REAL(8), intent(in) :: x
  REAL(8), intent(out) :: jl(0:lmax)
  ! local variables
  ! staring value for l above lmax (suitable for lmax < 50)
  INTEGER, parameter :: lst=25
  INTEGER l
  ! rescale limit
  REAL(8), parameter :: rsc=1.d150,rsci=1.d0/rsc
  REAL(8) :: xi,sx,cx
  REAL(8) :: j0,j1,jt,t1
  IF((lmax.lt.0).or.(lmax.gt.50)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(sbessel): lmax out of range : ",I8)') lmax
    WRITE(*,*)
    stop
  ENDIF 
  IF((x.lt.0.d0).or.(x.gt.1.d5)) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(sbessel): x out of range : ",G18.10)') x
    WRITE(*,*)
    stop
  ENDIF 
  ! treat x << 1
  IF(x.lt.1.d-8) THEN 
    jl(0)=1.d0
    t1=1.d0
    DO l=1,lmax
      t1=t1*x/(2*l+1)
      jl(l)=t1
    ENDDO 
    RETURN 
  ENDIF 
  xi=1.d0/x
  sx=sin(x)
  cx=cos(x)
  jl(0)=sx*xi
  IF(lmax.eq.0) RETURN 
  jl(1)=(jl(0)-cx)*xi
  IF(lmax.eq.1) RETURN 
  ! for x < lmax recurse down
  IF(x.lt.lmax) THEN 
  ! start from truly random numbers
    j1=0.6370354636449841609d0*rsci
    j0=0.3532702964695481204d0*rsci
    DO l=lmax+lst,lmax+1,-1
      jt=(2*l+1)*j1*xi-j0
      j0=j1
      j1=jt
  ! check for overflow
      IF(abs(j1).gt.rsc) THEN 
  ! rescale
        jt=jt*rsci
        j0=j0*rsci
        j1=j1*rsci
      ENDIF 
    ENDDO 
    DO l=lmax,2,-1
      jt=(2*l+1)*j1*xi-j0
      j0=j1
      j1=jt
  ! check for overflow
      IF(abs(j1).gt.rsc) THEN 
  ! rescale
        jt=jt*rsci
        j0=j0*rsci
        j1=j1*rsci
        jl(l+1:lmax)=jl(l+1:lmax)*rsci
      ENDIF 
      jl(l)=j0
    ENDDO 
    j0=3*j1*xi-j0
  ! rescaling constant
    t1=1.d0/((j0-x*j1)*cx+x*j0*sx)
    jl(2:)=t1*jl(2:)
    RETURN 
  ELSE 
  ! for large x recurse up
    j0=jl(0)
    j1=jl(1)
    DO l=2,lmax
      jt=(2*l-1)*j1*xi-j0
      j0=j1
      j1=jt
      jl(l)=j1
    ENDDO 
    RETURN 
  ENDIF 
END SUBROUTINE 
