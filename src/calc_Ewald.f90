SUBROUTINE calc_Ewald()

  USE m_constants, ONLY : PI
  USE m_atoms, ONLY : Natoms, &
                      atm2species, &
                      AtomicValences, &
                      AtomicCoords
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
  INTEGER :: ia, ja, isp, jsp
  REAL(8) :: invLatVecs(3,3)
  REAL(8), ALLOCATABLE :: tau(:,:)

  ALLOCATE( tau(3,Natoms) )
  CALL inv_m3x3(LatVecs,invLatVecs)
  tau(:,:) = matmul(invLatVecs,AtomicCoords)

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

  gcut = 2.d0
  ebsl = 1d-8

  tpi = 2.d0*pi
  on = volcry/(4.d0*pi)
  con2 = (4.d0*pi)/volcry
  glast2 = gcut*gcut
  gexp = -log(ebsl)
  eta = glast2/gexp

  cccc = sqrt(eta/pi)

  x = 0.d0
  totalcharge = 0.d0
  DO ia = 1,Natoms
    isp = atm2species(ia)
    x = x + AtomicValences(isp)**2
    totalcharge = totalcharge + AtomicValences(isp)
  ENDDO

  ewald = -cccc*x - 4.d0*pi*(totalcharge**2)/(volcry*eta)

  tmax = sqrt(2.d0*gexp/eta)
  seta = sqrt(eta)/2.d0

  mmm1 = nint(tmax/t1m + 1.5d0)
  mmm2 = nint(tmax/t2m + 1.5d0)
  mmm3 = nint(tmax/t3m + 1.5d0)
      
  DO ia = 1,Natoms
  DO ja = 1,Natoms
    v(:) = (tau(1,ia)-tau(1,ja))*t1(:) + (tau(2,ia)-tau(2,ja))*t2(:) &
         + (tau(3,ia)-tau(3,ja))*t3(:)
    isp = atm2species(ia)
    jsp = atm2species(ja)
    prd = AtomicValences(isp)*AtomicValences(jsp)
    DO i = -mmm1, mmm1
    DO j = -mmm2, mmm2
    DO k = -mmm3, mmm3
      IF( (ia /= ja) .or. ( (abs(i) + abs(j) + abs(k)) /= 0) ) THEN 
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

  mmm1 = nint(gcut/g1m + 1.5d0)
  mmm2 = nint(gcut/g2m + 1.5d0)
  mmm3 = nint(gcut/g3m + 1.5d0)
  
  DO i = -mmm1, mmm1
  DO j = -mmm2, mmm2
  DO k = -mmm3, mmm3
    IF( (abs(i)+abs(j)+abs(k)).ne.0 ) THEN 
      w(:) = i*g1(:) + j*g2(:) + k*g3(:)
      rmag2 = DOT_PRODUCT(w,w)
      x = con2*exp(-rmag2/eta)/rmag2
      DO ia = 1,Natoms
      DO ja = 1,Natoms
        v(:) = tau(:,ia) - tau(:,ja)
        isp = atm2species(ia)
        jsp = atm2species(ja)
        prd = AtomicValences(isp)*AtomicValences(jsp)
        arg = tpi*(i*v(1) + j*v(2) + k*v(3))
        ewald = ewald + x*prd*cos(arg)
      ENDDO 
      ENDDO 
    ENDIF 
  ENDDO 
  ENDDO 
  ENDDO

  E_nn = ewald/2.d0  ! convert to Ry to Ha

  DEALLOCATE( tau )

END SUBROUTINE 
