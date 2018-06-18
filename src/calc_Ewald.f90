SUBROUTINE calc_Ewald()

  USE m_constants, ONLY : PI
  USE m_atoms, ONLY : Nion => Natoms, &
                      ntyp => Nspecies, &
                      ityp => atm2species, &
                      q => AtomicValences, &
                      tau => AtomicCoords
  USE m_cell, ONLY : LatVecs
  USE m_energies, ONLY : E_nn

  IMPLICIT NONE 
  REAL(8) :: t1(3), t2(3), t3(3)
  REAL(8) :: g1(3), g2(3), g3(3)
  REAL(8) :: v(3), w(3)
  REAL(8) :: volcry, x, tpi, totalcharge, tmax, t1m, t2m, t3m
  REAL(8) :: seta, rmag2, prd, on, glast2, gexp, gcut, g1m, g2m, g3m
  REAL(8) :: cccc, ewald, eta, con2, ebsl, arg
  INTEGER :: mmm1, mmm2, mmm3, i, j, k
  INTEGER :: a, b, isp

  t1(:) = LatVecs(:,1)
  t2(:) = LatVecs(:,2)
  t3(:) = LatVecs(:,3)

  volcry   = t1(1)*(t2(2)*t3(3)-t2(3)*t3(2)) +  &
             t1(2)*(t2(3)*t3(1)-t2(1)*t3(3)) +  &
             t1(3)*(t2(1)*t3(2)-t2(2)*t3(1))

  g1(1) = 2 * pi * (t2(2)*t3(3)-t2(3)*t3(2))/volcry
  g1(2) = 2 * pi * (t2(3)*t3(1)-t2(1)*t3(3))/volcry
  g1(3) = 2 * pi * (t2(1)*t3(2)-t2(2)*t3(1))/volcry
  g2(1) = 2 * pi * (t3(2)*t1(3)-t3(3)*t1(2))/volcry
  g2(2) = 2 * pi * (t3(3)*t1(1)-t3(1)*t1(3))/volcry
  g2(3) = 2 * pi * (t3(1)*t1(2)-t3(2)*t1(1))/volcry
  g3(1) = 2 * pi * (t1(2)*t2(3)-t1(3)*t2(2))/volcry
  g3(2) = 2 * pi * (t1(3)*t2(1)-t1(1)*t2(3))/volcry
  g3(3) = 2 * pi * (t1(1)*t2(2)-t1(2)*t2(1))/volcry
  
  volcry = abs(volcry)

  t1m = SQRT(DOT_PRODUCT(t1,t1))
  t2m = SQRT(DOT_PRODUCT(t2,t2))
  t3m = SQRT(DOT_PRODUCT(t3,t3))
  g1m = SQRT(DOT_PRODUCT(g1,g1))
  g2m = SQRT(DOT_PRODUCT(g2,g2))
  g3m = SQRT(DOT_PRODUCT(g3,g3))

  !ALLOCATE( q(nion), tau(3,nion) )

!  Write(6,*) 'For each ion, input q and tau (in fractional coordinates of T)'

  !DO i=1,nion
  !  READ(5,*) q(i),tau(1,i),tau(2,i),tau(3,i)
  !ENDDO

!  Write(6,*) ' Input Gcut for reciprocal lattice sum and error tolerance'
!  Read(5,*) gcut , ebsl
  gcut = 2.d0
  ebsl = 1d-8

  tpi = 2.d0*pi
  on = volcry/(4.d0*pi)
  con2 = (4.d0*pi)/volcry
  glast2 = gcut*gcut
  gexp = -log(ebsl)
  eta = glast2/gexp

  WRITE(*,*) 'eta value for this calculation' , eta
  
  cccc = sqrt(eta/pi)

  x = 0.d0
  totalcharge = 0.d0
  DO i = 1,nion 
    isp = ityp(i)
    x = x + q(isp)**2
    totalcharge = totalcharge + q(isp)
  ENDDO

  totalcharge = sum(q) 
  
  WRITE(*,*) 'Total charge = ', totalcharge
  
  ewald = -cccc*x - 4.d0*pi*(totalcharge**2)/(volcry*eta)

  tmax = sqrt(2.d0*gexp/eta)
  seta = sqrt(eta)/2.d0

  mmm1 = tmax/t1m + 1.5d0
  mmm2 = tmax/t2m + 1.5d0
  mmm3 = tmax/t3m + 1.5d0  
      
  WRITE(*,*) 'Lattice summation indices -- ', mmm1,mmm2,mmm3
  DO a = 1,Nion
  DO b = 1,Nion
    v(:) = (tau(1,a)-tau(1,b))*t1(:) + (tau(2,a)-tau(2,b))*t2(:) &
         + (tau(3,a)-tau(3,b))*t3(:)
    prd = q(a)*q(b)
    DO i = -mmm1, mmm1
    DO j = -mmm2, mmm2
    DO k = -mmm3, mmm3
      IF( (a.ne.b).or.((abs(i)+abs(j)+abs(k)).ne.0) ) THEN 
        w(:) = v(:) + i*t1 + j*t2 + k*t3
        rmag2 = sqrt(DOT_PRODUCT(w,w))
        arg = rmag2*seta 
        ewald = ewald + prd*erfc(arg)/rmag2
      ENDIF 
    ENDDO 
    ENDDO 
    ENDDO 
  ENDDO 
  ENDDO 

  mmm1 = gcut/g1m + 1.5d0
  mmm2 = gcut/g2m + 1.5d0
  mmm3 = gcut/g3m + 1.5d0
      
  WRITE(*,*) 'Reciprocal lattice summation indices --', mmm1,mmm2,mmm3
  DO i = -mmm1, mmm1
  DO j = -mmm2, mmm2
  DO k = -mmm3, mmm3
    IF( (abs(i)+abs(j)+abs(k)).ne.0 ) THEN 
      w(:) = i*g1(:) + j*g2(:) + k*g3(:)
      rmag2 = DOT_PRODUCT(w,w)
      x = con2*exp(-rmag2/eta)/rmag2
      DO a = 1,nion
      DO b = 1,nion
        v(:) = tau(:,a) - tau(:,b)
        prd = q(a)*q(b)
        arg = tpi*(i*v(1) + j*v(2) + k*v(3))
        ewald = ewald + x*prd*cos(arg)
      ENDDO 
      ENDDO 
    ENDIF 
  ENDDO 
  ENDDO 
  ENDDO

  WRITE(*,'(1x,A,F18.10)') 'Ewald energy in Ry', ewald
  WRITE(*,'(1x,A,F18.10)') 'Ewald energy in Ha', ewald/2.d0

  E_nn = ewald/2.d0

END SUBROUTINE 
