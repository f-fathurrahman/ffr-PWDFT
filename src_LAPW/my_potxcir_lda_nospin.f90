SUBROUTINE my_potxcir_lda_nospin(xctype_, &
  rhoir_, &
  exir_, ecir_, &
  vxcir_)
  !
  USE m_gvectors, ONLY: ngtot
  USE m_spin, ONLY: spinpol
  USE m_density_pot_xc, ONLY: xcgrad
  USE m_states, ONLY: swidth
  USE m_oep_hf, ONLY: hybridc, hybrid
  USE modxcifc, ONLY: xcifc
  !
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: xctype_(3)
  REAL(8), intent(in) :: rhoir_(ngtot)
  REAL(8), intent(out) :: exir_(ngtot),ecir_(ngtot)
  REAL(8), intent(out) :: vxcir_(ngtot)
  ! local variables
  INTEGER :: n
  REAL(8) :: t1
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: vx(:)
  REAL(8), ALLOCATABLE :: vc(:)
  
  if( xcgrad /= 0 ) then
    stop 'my_potxcir_lda_nospin: xcgrad /= 0 is not supported in this subroutine'
  endif

  IF(spinpol) THEN 
    stop 'my_potxcir_lda_nospin: spinpol == true is not supported in this subroutine'
  ENDIF 

  n = ngtot
  ALLOCATE(vx(n), vc(n))

  !--------------------------!
  !     spin-unpolarised     !
  !--------------------------!
  CALL xcifc(xctype_, n=n, tempa=swidth, rho=rhoir_, ex=exir_, ec=ecir_, vx=vx, vc=vc)

  ! hybrid functionals
  IF(hybrid) THEN 
    t1 = 1.d0 - hybridc
    ! scale exchange part of energy
    exir_(:) = t1*exir_(:)
    ! scale exchange part of potential
    vxcir_(:) = t1*vx(:) + vc(:)
  ELSE 
    vxcir_(:) = vx(:) + vc(:)
  ENDIF 


  DEALLOCATE(vx,vc)

  RETURN 
END SUBROUTINE 

