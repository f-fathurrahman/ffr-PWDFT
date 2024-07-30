! !INPUT/OUTPUT PARAMETERS:
!   sol   : speed of light in atomic units (in,real)
!   l     : angular momentum quantum number (in,integer)
!   nr    : number of radial mesh points (in,integer)
!   r     : radial mesh (in,real(nr))
!   vr    : potential on radial mesh (in,real(nr))
!   eps   : energy search tolerance (in,real)
!   demax : maximum allowed change from the input energy; enforced only if e < 0
!           (in,real)
!   e     : input energy and returned band energy (inout,real)
!   fnd   : set to .true. if the band energy is found (out,logical)

!------------------------------------------------------------
SUBROUTINE my_findband(sol, l, nr, r, vr, eps, demax, e, fnd)
!------------------------------------------------------------
  IMPLICIT NONE
  ! arguments
  REAL(8), INTENT(in) :: sol
  INTEGER, INTENT(in) :: l,nr
  REAL(8), INTENT(in) :: r(nr),vr(nr)
  REAL(8), INTENT(in) :: eps,demax
  REAL(8), INTENT(inout) :: e
  LOGICAL, INTENT(out) :: fnd
  ! local variables
  LOGICAL :: ft,fb
  ! maximum number of steps
  INTEGER, PARAMETER :: maxstp=250
  INTEGER :: ip,ie,nn
  ! initial step size
  REAL(8), PARAMETER :: de0=0.001d0
  REAL(8) :: de,et,eb,t,tp
  ! automatic arrays
  REAL(8) :: p0(nr),p1(nr),q0(nr),q1(nr)
  !
  ft = .false.
  fb = .false.
  fnd = .false.
  et = e
  eb = e
  ! two-pass loop
  DO ip = 1,2
    ! find the top of the band
    tp = 0.d0
    de = de0
    do ie = 1,maxstp
      et = et + de
      IF( e < 0.d0 ) THEN 
        IF( et > e+demax ) EXIT 
      ENDIF
      CALL my_rschrodint(sol,l,et,nr,r,vr,nn,p0,p1,q0,q1)
      t = p0(nr)
      IF( ie > 1 ) THEN
        IF( t*tp <= 0.d0 ) THEN 
          IF( abs(de) < eps ) THEN 
            IF(fb) GOTO 10
            ft = .true.
            eb = et + 5.d0*abs(de0)
            EXIT
          ENDIF
          de = -0.5d0*de
        ELSE 
          de = 1.5d0*de
        ENDIF
      ENDIF
      tp = t
    ENDDO
    IF( fb ) return
    ! find the bottom of the band
    tp = 0.d0
    de = -de0
    DO ie = 1,maxstp
      eb = eb + de
      IF( eb < e-demax ) RETURN 
      CALL my_rschrodint(sol, l, eb, nr, r, vr, nn, p0, p1, q0, q1)
      t = p1(nr)
      IF( ie > 1) THEN 
        IF( t*tp <= 0.d0) THEN 
          IF( abs(de) < eps) THEN 
            IF( ft ) GOTO 10
            fb = .true.
            et = eb - 5.d0*abs(de0)
            EXIT 
          ENDIF
          de = -0.5d0*de
        ELSE 
          de = 1.5d0*de
        ENDIF
      ENDIF
      tp = t
    ENDDO
  ENDDO
  RETURN 
  !
  10 CONTINUE
  ! set the band energy halfway between top and bottom
  e = (et + eb)/2.d0
  fnd = .true.
  RETURN 
END SUBROUTINE 


