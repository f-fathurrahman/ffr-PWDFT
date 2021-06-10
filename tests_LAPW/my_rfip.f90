SUBROUTINE my_rfip(ip, np, vpl, zfft)
  USE m_atoms, ONLY: natoms, natmtot, idxas, rsp, nspecies, atposc
  USE m_gvectors, ONLY: ngtot, vgc, igfft, ngvec
  USE m_lattice, ONLY: epslat, avec
  USE m_muffin_tins, ONLY: lmmaxo, nrmt, nrmti, nrmtmax, rmt, lmaxo, lmaxi
  IMPLICIT NONE 
  ! arguments
  INTEGER :: ip
  INTEGER :: np
  REAL(8) :: vpl(3,np)
  COMPLEX(8) :: zfft(ngtot)
  REAL(8) :: rfmt1(lmmaxo,nrmtmax,natmtot)
  REAL(8) :: fp(np)
  ! local variables
  INTEGER :: is,ia,ias,nr,nri
  INTEGER :: ir0,ir,lmax,l,m,lm
  INTEGER :: ig,ifg,i1,i2,i3,i,j
  REAL(8) :: rmt2,r,sum,ya(4),t1
  REAL(8) :: v1(3),v2(3),v3(3),v4(3),v5(3)
  ! automatic arrays
  REAL(8) rlm(lmmaxo)
  
  v2(:) = vpl(:,ip)
  CALL r3frac(epslat,v2)
  ! convert point to Cartesian coordinates
  CALL r3mv(avec,v2,v1)

  ! check if point is in a muffin-tin
  DO is=1,nspecies
    nr=nrmt(is)
    nri=nrmti(is)
    rmt2=rmt(is)**2
    DO ia=1,natoms(is)
      ias=idxas(ia,is)
      v2(:)=v1(:)-atposc(:,ia,is)
      DO i1=-1,1
        v3(:)=v2(:)+dble(i1)*avec(:,1)
        DO i2=-1,1
          v4(:)=v3(:)+dble(i2)*avec(:,2)
          DO i3=-1,1
            v5(:)=v4(:)+dble(i3)*avec(:,3)
            t1=v5(1)**2+v5(2)**2+v5(3)**2
            !
            IF(t1 < rmt2) THEN 
              r=sqrt(t1)
              CALL genrlmv(lmaxo,v5,rlm)
              !
              DO ir=1,nr
                !
                IF(rsp(ir,is) >= r) THEN 
                  write(*,*)
                  write(*,*) 'inside the muffin tin'
                  write(*,*) 'r = ', r
                  write(*,*) 
                  !
                  IF(ir.le.3) THEN 
                    ir0=1
                  ELSEIF(ir.gt.nr-2) THEN 
                    ir0=nr-3
                  else
                    ir0=ir-2
                  ENDIF 
                  !
                  r=max(r,rsp(1,is))
                  IF(ir0 <= nri) THEN 
                    lmax=lmaxi
                  else
                    lmax=lmaxo
                  ENDIF 
                  !
                  sum=0.d0
                  lm=0
                  DO l=0,lmax
                    DO m=-l,l
                      lm=lm+1
                      DO j=1,4
                        i=ir0+j-1
                        ya(j)=rfmt1(lm,i,ias)
                      ENDDO 
                      t1=poly4(rsp(ir0,is),ya,r)
                      sum=sum+t1*rlm(lm)
                    ENDDO 
                  ENDDO 
                  goto 10 ! move out, set fp, return
                ENDIF  ! rsp(ir,is) >= r
              ENDDO  ! loop over all ir in radial grid of current atom
            ENDIF  ! if t1 < rmt2
          ENDDO 
        ENDDO 
      ENDDO 
    ENDDO  ! loop over atoms in species
  ENDDO  ! loop over nspecies
  !
  ! otherwise use direct Fourier transform of interstitial function
  !
  sum=0.d0
  DO ig=1,ngvec
    ifg = igfft(ig)
    t1 = vgc(1,ig)*v1(1) + vgc(2,ig)*v1(2) + vgc(3,ig)*v1(3)
    sum = sum + dble(zfft(ifg)*cmplx(cos(t1),sin(t1),8))
  ENDDO 
  10 continue
  fp(ip)=sum
  RETURN 


CONTAINS

!------------------------------
REAL(8) function poly4(xa,ya,x)
!------------------------------
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(in) :: xa(4),ya(4),x
  ! local variables
  REAL(8) x0,x1,x2,x3,y0,y1,y2,y3
  REAL(8) c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
  ! evaluate the polynomial coefficients
  x0=xa(1)
  x1=xa(2)-x0; x2=xa(3)-x0; x3=xa(4)-x0
  t4=x1-x2; t5=x1-x3; t6=x2-x3
  y0=ya(1)
  y1=ya(2)-y0; y2=ya(3)-y0; y3=ya(4)-y0
  t1=x1*x2*y3; t2=x2*x3*y1; t3=x1*x3
  t0=1.d0/(x2*t3*t4*t5*t6)
  t3=t3*y2
  c3=t1*t4+t2*t6-t3*t5
  t4=x1**2; t5=x2**2; t6=x3**2
  c2=t1*(t5-t4)+t2*(t6-t5)+t3*(t4-t6)
  c1=t1*(x2*t4-x1*t5)+t2*(x3*t5-x2*t6)+t3*(x1*t6-x3*t4)
  t1=x-x0
  ! evaluate the polynomial
  poly4=y0+t0*t1*(c1+t1*(c2+c3*t1))
  RETURN 
END FUNCTION

END SUBROUTINE 