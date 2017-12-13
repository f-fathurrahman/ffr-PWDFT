! from QE-6.0

SUBROUTINE gen_ylmr2(lmax2, NR, R, Rm, ylm)
  !
  !     Real spherical harmonics ylm(R) up to l=lmax
  !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
  !     Numerical recursive algorithm based on the one given in Numerical
  !     Recipes but avoiding the calculation of factorials that generate
  !     overflow for lmax > 11
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = 8
  REAL(DP), PARAMETER :: PI = 4.d0*atan(1.d0)
  REAL(DP), PARAMETER :: FPI = 4.d0*PI
  !
  INTEGER, INTENT(in) :: lmax2, NR
  REAL(DP), INTENT(in) :: R(3, NR), Rm(NR)
  !
  ! BEWARE: Rm = R(1)^2 + R(2)^2 + R(3)^2  is not checked on input
  !         incorrect results will ensue if the above does not hold
  !
  REAL(DP), INTENT(out) :: ylm(NR,lmax2)
  !
  ! local variables
  !
  REAL(DP), PARAMETER :: eps = 1.0d-9
  REAL(DP), ALLOCATABLE :: cost(:), sent(:), phi(:), Q(:,:,:)
  REAL(DP) :: c, gmod
  INTEGER :: lmax, ir, l, m, lm
  !

  IF( NR < 1 .or. lmax2 < 1 ) RETURN

  DO lmax = 0, 25
    WRITE(*,*) 'lmax = ', lmax
    IF ((lmax+1)**2 == lmax2) GOTO 10
  ENDDO
  WRITE(*,*) 'Error ylmr: l > 25 or wrong number of Ylm required, lmax, lmax2 = ', lmax, lmax2
  STOP

10 CONTINUE
  !
  IF(lmax == 0) THEN
    ylm(:,1) =  sqrt(1.d0 / fpi)
    RETURN
  ENDIF

  !
  !  theta and phi are polar angles, cost = cos(theta)
  !
  ALLOCATE( cost(NR), sent(NR), phi(NR), Q(NR,0:lmax,0:lmax) )
  Q(:,:,:) = 0.d0
  !

  DO ir = 1, NR
     gmod = sqrt(Rm(ir))
     IF( gmod < eps ) THEN
        cost(ir) = 0.d0
     ELSE
        cost(ir) = R(3,ir)/gmod
     ENDIF 
     !
     !  beware the arc tan, it is defined modulo pi
     !
     IF( R(1,ir) > eps ) THEN
        phi (ir) = atan( R(2,ir)/R(1,ir) )
     ELSEIF( R(1,ir) < -eps) THEN
        phi(ir) = atan( R(2,ir)/R(1,ir) ) + PI
     ELSE
        phi (ir) = sign( PI/2.d0,R(2,ir) )
     ENDIF
     sent(ir) = sqrt(max(0d0,1.d0-cost(ir)**2))
     !WRITE(*,'(1x,I5,3F18.10)') ir, phi(ir), cost(ir), sent(ir)
  ENDDO
  !
  !  Q(:,l,m) are defined as sqrt ((l-m)!/(l+m)!) * P(:,l,m) where
  !  P(:,l,m) are the Legendre Polynomials (0 <= m <= l)
  !
  lm = 0
  DO l = 0, lmax

    !WRITE(*,*) 'l = ', l
    c = sqrt (DBLE(2*l+1) / fpi)
    
    IF( l == 0 ) THEN
      DO ir = 1, NR
        Q(ir,0,0) = 1.d0
      ENDDO

    ELSEIF ( l == 1 ) THEN
      DO ir = 1, NR
        Q(ir,1,0) =  cost(ir)
        Q(ir,1,1) = -sent(ir)/sqrt(2.d0)
      ENDDO 

    ELSE
      !
      !  recursion on l for Q(:,l,m)
      !
      DO m = 0, l - 2
        !write(*,*) 'Enter this ....'
        DO ir = 1, NR
          Q(ir,l,m) = cost(ir)*(2*l-1)/sqrt(DBLE(l*l-m*m))*Q(ir,l-1,m) &
                       - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*Q(ir,l-2,m)
        ENDDO
        !write(*,*) 'sum(Q) = ', sum(Q(:,l,m))
      ENDDO
      
      DO ir = 1, NR
        Q(ir,l,l-1) = cost(ir) * sqrt(DBLE(2*l-1)) * Q(ir,l-1,l-1)
      ENDDO
      
      DO ir = 1, NR
        Q(ir,l,l)   = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent(ir)*Q(ir,l-1,l-1)
      ENDDO

    ENDIF 
    !
    ! Y_lm, m = 0
    !
    lm = lm + 1
    !WRITE(*,*) 'm = 0: lm = ', lm
    DO ir = 1, NR
      ylm(ir, lm) = c * Q(ir,l,0)
    ENDDO
    !
    DO m = 1, l
      !
      ! Y_lm, m > 0
      !
      lm = lm + 1
      !WRITE(*,*) 'm > 0: lm = ', lm
      DO ir = 1, NR
        ylm(ir, lm) = c * sqrt(2.d0) * Q(ir,l,m) * cos (m*phi(ir))
      ENDDO
      !
      ! Y_lm, m < 0
      !
      lm = lm + 1
      !WRITE(*,*) 'm < 0: lm = ', lm
      DO ir = 1, NR
        ylm(ir, lm) = c * sqrt(2.d0) * Q(ir,l,m) * sin (m*phi(ir))
      ENDDO

    ENDDO

  ENDDO

  !
  DEALLOCATE(cost, sent, phi, Q)
  !
  RETURN
END SUBROUTINE 
