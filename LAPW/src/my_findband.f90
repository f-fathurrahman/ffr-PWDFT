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
subroutine my_findband(sol, l, nr, r, vr, eps, demax, e, fnd)
!------------------------------------------------------------
  implicit none
  ! arguments
  real(8), intent(in) :: sol
  integer, intent(in) :: l,nr
  real(8), intent(in) :: r(nr),vr(nr)
  real(8), intent(in) :: eps,demax
  real(8), intent(inout) :: e
  logical, intent(out) :: fnd
  ! local variables
  logical ft,fb
  ! maximum number of steps
  integer, parameter :: maxstp=250
  integer ip,ie,nn
  ! initial step size
  real(8), parameter :: de0=0.001d0
  real(8) de,et,eb,t,tp
  ! automatic arrays
  real(8) p0(nr),p1(nr),q0(nr),q1(nr)
  ft=.false.
  fb=.false.
  fnd=.false.
  et=e
  eb=e
  ! two-pass loop
  do ip=1,2
  ! find the top of the band
    tp=0.d0
    de=de0
    do ie=1,maxstp
      et=et+de
      if (e.lt.0.d0) then
        if (et.gt.e+demax) exit
      ENDIF
      call rschrodint(sol,l,et,nr,r,vr,nn,p0,p1,q0,q1)
      t=p0(nr)
      if (ie.gt.1) then
        if (t*tp.le.0.d0) then
          if (abs(de) .lt. eps) then
            if(fb) goto 10
            ft = .true.
            eb = et + 5.d0*abs(de0)
            exit
          ENDIF
          de=-0.5d0*de
        else
          de=1.5d0*de
        ENDIF
      ENDIF
      tp=t
    ENDDO
    if (fb) return
  ! find the bottom of the band
    tp=0.d0
    de=-de0
    do ie=1,maxstp
      eb=eb+de
      if (eb.lt.e-demax) return
      call my_rschrodint(sol,l,eb,nr,r,vr,nn,p0,p1,q0,q1)
      t=p1(nr)
      if (ie.gt.1) then
        if (t*tp.le.0.d0) then
          if (abs(de).lt.eps) then
            if (ft) goto 10
            fb=.true.
            et=eb-5.d0*abs(de0)
            exit
          ENDIF
          de=-0.5d0*de
        else
          de=1.5d0*de
        ENDIF
      ENDIF
      tp=t
    ENDDO
  ENDDO
  return
  !
  10 continue
  ! set the band energy halfway between top and bottom
  e = (et + eb)/2.d0
  fnd = .true.
  return
end subroutine
