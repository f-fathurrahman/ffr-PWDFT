! !INPUT/OUTPUT PARAMETERS:
!   nr     : number of radial mesh points (in,integer)
!   nri    : number of points on inner part of muffin-tin (in,integer)
!   ld     : leading dimension (in,integer)
!   rl     : r^l on the radial mesh (in,real(ld,-lmaxo-1:lmaxo+2))
!   wpr    : weights for partial integration on radial mesh (in,real(4,nr))
!   zrhomt : muffin-tin charge density (in,complex(*))
!   zvclmt : muffin-tin Coulomb potential (out,complex(*))
SUBROUTINE my_zpotclmt(nr,nri,ld,rl,wpr,zrhomt,zvclmt)
  USE m_constants, ONLY: fourpi
  USE m_muffin_tins, ONLY: lmaxo, lmmaxo, lmmaxi, lmaxi
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr,nri,ld
  REAL(8), intent(in) :: rl(ld,-lmaxo-1:lmaxo+2),wpr(4,nr)
  COMPLEX(8), intent(in) :: zrhomt(*)
  COMPLEX(8), intent(out) :: zvclmt(*)
  ! local variables
  INTEGER :: nro,iro,ir
  INTEGER :: l,l1,l2,l3
  INTEGER :: m,lm,npi,i
  REAL(8) :: r1,r2,t0,t1,t2,t3,t4
  ! automatic arrays
  REAL(8) :: f1(nr),f2(nr),f3(nr),f4(nr),f5(nr)
  complex(8) :: sz

  nro = nr - nri
  iro = nri + 1
  npi = lmmaxi*nri
  lm = 0
  write(*,*) 'lmaxi = ', lmaxi
  write(*,*) 'my_zpotclmt: sum(wpr) = ', sum(wpr)
  write(*,*) 'my_zpotclmt: sum(rl) = ', sum(rl)

  sz = cmplx(0.d0,0.d0,kind=8)
  DO l = 0,lmaxi
    l1 =  l + 2
    l2 = -l + 1
    l3 = -l - 1
    write(*,'(1x,A,3I4)') 'l1 l2 l3 = ', l1, l2, l3
    t0 = fourpi/dble(2*l+1)
    DO m = -l,l
      lm = lm + 1
      i = lm
      DO ir = 1,nri
        t1 = dble(zrhomt(i))
        t2 = aimag(zrhomt(i))
        r1 = rl(ir,l1)
        r2 = rl(ir,l2)
        f1(ir) = t1*r1
        f2(ir) = t2*r1
        f3(ir) = t1*r2
        f4(ir) = t2*r2
        !write(*,'(1x,A,I8,2F18.10)') 'zrhomt: ', i, t1, t2
        i = i + lmmaxi
      ENDDO 
      DO ir = iro,nr
        t1 = dble(zrhomt(i))
        t2 = aimag(zrhomt(i))
        r1 = rl(ir,l1)
        r2 = rl(ir,l2)
        f1(ir) = t1*r1
        f2(ir) = t2*r1
        f3(ir) = t1*r2
        f4(ir) = t2*r2
        !write(*,'(1x,A,I8,2F18.10)') 'zrhomt: ', i, t1, t2
        i = i + lmmaxo
      ENDDO 
      CALL splintwp(nr,wpr,f1,f5)
      CALL splintwp(nr,wpr,f2,f1)
      CALL splintwp(nr,wpr,f3,f2)
      CALL splintwp(nr,wpr,f4,f3)
      t1 = f2(nr)
      t2 = f3(nr)
      i = lm
      DO ir = 1,nri
        r1 = t0*rl(ir,l3)
        r2 = t0*rl(ir,l)
        t3 = r1*f5(ir) + r2*(t1 - f2(ir))
        t4 = r1*f1(ir) + r2*(t2 - f3(ir))
        zvclmt(i) = cmplx(t3,t4,8)
        !write(*,'(1x,A,I8,2F18.10)') 'zvclmt: ', i, t3, t4
        sz = sz + zvclmt(i)
        i = i + lmmaxi
      ENDDO 
      DO ir = iro,nr
        r1 = t0*rl(ir,l3)
        r2 = t0*rl(ir,l)
        t3 = r1*f5(ir) + r2*(t1 - f2(ir))
        t4 = r1*f1(ir) + r2*(t2 - f3(ir))
        zvclmt(i) = cmplx(t3,t4,8)
        !write(*,'(1x,A,I8,2F18.10)') 'zvclmt: ', i, t3, t4
        sz = sz + zvclmt(i)
        i = i + lmmaxo
      ENDDO 
    ENDDO 
  ENDDO 
  write(*,*) 'After lmaxi sz = ', sz
  !stop 'ffr 99'


  write(*,*) 'lmaxo = ', lmaxo
  DO l=lmaxi+1,lmaxo
    l1 =  l + 2
    l2 = -l + 1
    l3 = -l - 1
    write(*,'(1x,A,3I4)') 'l1 l2 l3 = ', l1, l2, l3
    t0=fourpi/dble(2*l+1)
    DO m=-l,l
      lm=lm+1
      i=npi+lm
      DO ir=iro,nr
        t1=dble(zrhomt(i))
        t2=aimag(zrhomt(i))
        r1=rl(ir,l1)
        r2=rl(ir,l2)
        f1(ir)=t1*r1
        f2(ir)=t2*r1
        f3(ir)=t1*r2
        f4(ir)=t2*r2
        i=i+lmmaxo
      ENDDO 
      CALL splintwp(nro,wpr(1,iro),f1(iro),f5(iro))
      CALL splintwp(nro,wpr(1,iro),f2(iro),f1(iro))
      CALL splintwp(nro,wpr(1,iro),f3(iro),f2(iro))
      CALL splintwp(nro,wpr(1,iro),f4(iro),f3(iro))
      t1=f2(nr); t2=f3(nr)
      i=npi+lm
      DO ir=iro,nr
        r1 = t0*rl(ir,l3); r2=t0*rl(ir,l)
        t3 = r1*f5(ir)+r2*(t1-f2(ir))
        t4 = r1*f1(ir)+r2*(t2-f3(ir))
        zvclmt(i) = cmplx(t3,t4,8)
        i = i + lmmaxo
      ENDDO 
    ENDDO 
  ENDDO 

  RETURN 
END SUBROUTINE 