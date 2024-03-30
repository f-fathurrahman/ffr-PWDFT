PROGRAM main
  IMPLICIT NONE 

  ! read input and initialize some variables
  CALL read_input()
  CALL init0()
  CALL init1()
  
  CALL info_crystal()
  CALL info_symmetry()
  CALL writesym()
  CALL info_gvectors()
  CALL info_muffin_tins()

  CALL rhoinit()

  CALL test_plot_around0()

END PROGRAM

!-----------------------------
subroutine test_plot_around0()
!-----------------------------
  USE m_density_pot_xc, ONLY: rhomt
  USE m_atoms, ONLY: natoms, natmtot, idxis, idxas, nspecies, atposc
  USE m_muffin_tins, ONLY: lmmaxo, nrmt, nrmti, nrmtmax, rmt, lmaxi, lmaxo
  USE m_lattice, ONLY: epslat, avec
  USE m_atomic_species, ONLY: rsp
  !
  implicit none
  !
  REAL(8), ALLOCATABLE :: rhomt_unpacked(:,:,:)
  INTEGER :: is,ia,ias,nr,nri
  INTEGER :: ir0,ir,lmax,l,m,lm
  INTEGER :: i1,i2,i3,i,j
  REAL(8) :: rmt2, r, sum, ya(4), t1
  REAL(8) :: v1(3), v2(3), v3(3), v4(3), v5(3)
  ! automatic arrays
  REAL(8) :: rlm(lmmaxo)


  ! First, we unpack the muffin-tin function
  ALLOCATE( rhomt_unpacked(lmmaxo,nrmtmax,natmtot) )
  DO ias = 1,natmtot
    is = idxis(ias)
    CALL rfmtpack(.false., nrmt(is), nrmti(is), rhomt(:,ias), rhomt_unpacked(:,:,ias))
  ENDDO 

  !
  ! adapted from rfip
  !

  !v2(:) = (/ 0.000001d0, 0.d0, 0.d0 /)
  !v2(:) = atposc(:,1,1) + (/ 0.0001d0, 0.d0, 0.d0 /)
  CALL r3frac(epslat, v2)
  CALL r3mv(avec, v2, v1)   ! convert point to Cartesian coordinates


  v1(:) = atposc(:,1,1) + (/ 0.0001d0, 0.d0, 0.d0 /)
  
  write(*,*) 'v1 = ', v1

  ! check if point is in a muffin-tin
  DO is = 1,nspecies
    nr = nrmt(is)
    nri = nrmti(is)
    rmt2 = rmt(is)**2
    DO ia = 1,natoms(is)
      write(*,*) 'atposc = ', atposc(:,ia,is)
      ias = idxas(ia,is)
      v2(:) = v1(:) - atposc(:,ia,is)
      ! we will loop over nearest neighboring periodic cells
      DO i1 = -1,1
        v3(:) = v2(:) + dble(i1)*avec(:,1)
        DO i2 = -1,1
          v4(:) = v3(:) + dble(i2)*avec(:,2)
          DO i3 = -1,1
            v5(:) = v4(:) + dble(i3)*avec(:,3)
            !
            t1 = v5(1)**2 + v5(2)**2 + v5(3)**2
            write(*,'(1x,3I5,F18.10)') i1, i2, i3, sqrt(t1)
            !
            IF(t1 < rmt2) THEN 
              r = sqrt(t1)
              CALL genrlmv(lmaxo,v5,rlm) ! generate real spherical harmonics
              !
              DO ir = 1,nr
                !
                IF(rsp(ir,is) >= r) THEN 
                  !write(*,*) 'inside the muffin tin'
                  !write(*,'(1x,A,ES18.10)') 'r = ', r
                  !
                  IF(ir .le. 3) THEN 
                    ir0 = 1
                  ELSEIF(ir .gt. nr-2) THEN 
                    ir0 = nr - 3
                  ELSE 
                    ir0 = ir - 2
                  ENDIF 
                  !write(*,*) 'ir0 = ', ir0
                  !
                  r = max(r, rsp(1,is))
                  !write(*,'(1x,A,ES18.10)') 'r after max = ', r
                  IF(ir0 <= nri) THEN 
                    write(*,*) 'inner muffin tin'
                    lmax = lmaxi
                  ELSE 
                    write(*,*) 'outer muffin tin'
                    lmax = lmaxo
                  ENDIF 
                  write(*,*) 'lmax = ', lmax
                  !
                  sum = 0.d0
                  lm = 0
                  ! This is sum over Ylm
                  ! if the point falls into inner MT, lmax is lmaxi
                  ! if the point falls into outer MT, lmax is lmaxo
                  DO l = 0,lmax
                    DO m = -l,l
                      lm = lm + 1
                      DO j = 1,4
                        i = ir0 + j - 1
                        ya(j) = rhomt_unpacked(lm,i,ias)
                        !write(*,'(1x,2I4,ES18.10)') lm, i, ya(j)
                      ENDDO 
                      !write(*,*)
                      t1 = poly4(rsp(ir0,is), ya, r)
                      sum = sum + t1*rlm(lm)
                    ENDDO 
                  ENDDO 
                  write(*,*) 'sum = ', sum
                  goto 10 ! move out, set fp, return
                ENDIF  ! rsp(ir,is) >= r
              ENDDO  ! loop over all ir in radial grid of current atom
            ENDIF  ! if t1 < rmt2
          ENDDO 
        ENDDO 
      ENDDO 
    ENDDO  ! loop over atoms in species
  ENDDO  ! loop over nspecies


  write(*,*) 'Point is outside muffin tin'


10 continue
  write(*,*) 'Finished'
  deallocate(rhomt_unpacked)

contains

!-----------------------------
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

end subroutine

include 'my_plot3d.f90'
include 'my_plotpt3d.f90'
