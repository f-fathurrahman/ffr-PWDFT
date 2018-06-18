SUBROUTINE init_gvec()
  USE m_constants, ONLY : PI
  USE m_realspace, ONLY : Ns
  USE m_cell, ONLY : RecVecs
  USE m_PWGrid, ONLY : Gv, Gv2, Ng, ecutrho, idx_g2r
  IMPLICIT NONE
  !
  INTEGER :: i, j, k, ig, ii, jj, kk, ip
  REAL(8) :: Gvv(3)
  REAL(8) :: Gvv2
  ! Function
  INTEGER :: mm_to_nn

  Ng = calc_Ng()
  ALLOCATE( Gv(3,Ng) )
  ALLOCATE( Gv2(Ng) )
  ALLOCATE( idx_g2r(Ng) )

  ig = 0
  ip = 0
  DO k = 0, Ns(3)-1
  DO j = 0, Ns(2)-1
  DO i = 0, Ns(1)-1
    !
    ip = ip + 1
    !
    ii = mm_to_nn( i, Ns(1) )
    jj = mm_to_nn( j, Ns(2) )
    kk = mm_to_nn( k, Ns(3) )
    Gvv(1) = RecVecs(1,1)*ii + RecVecs(1,2)*jj + RecVecs(1,3)*kk
    Gvv(2) = RecVecs(2,1)*ii + RecVecs(2,2)*jj + RecVecs(2,3)*kk
    Gvv(3) = RecVecs(3,1)*ii + RecVecs(3,2)*jj + RecVecs(3,3)*kk
    !
    Gvv2 = Gvv(1)**2 + Gvv(2)**2 + Gvv(3)**2
    !
    IF( Gvv2 <= ecutrho ) THEN 
      ig = ig + 1
      Gv(:,ig) = Gvv(:)
      Gv2(ig) = Gvv2
      idx_g2r(ig) = ip
    ENDIF 
  ENDDO
  ENDDO
  ENDDO

CONTAINS 


FUNCTION calc_Ng() RESULT(Ng)
  IMPLICIT NONE 
  INTEGER :: Ng
  INTEGER :: i, j, k

  Ng = 0
  DO k = 0, Ns(3)-1
  DO j = 0, Ns(2)-1
  DO i = 0, Ns(1)-1
    !
    ii = mm_to_nn( i, Ns(1) )
    jj = mm_to_nn( j, Ns(2) )
    kk = mm_to_nn( k, Ns(3) )
    !
    Gvv(1) = RecVecs(1,1)*ii + RecVecs(1,2)*jj + RecVecs(1,3)*kk
    Gvv(2) = RecVecs(2,1)*ii + RecVecs(2,2)*jj + RecVecs(2,3)*kk
    Gvv(3) = RecVecs(3,1)*ii + RecVecs(3,2)*jj + RecVecs(3,3)*kk
    !
    Gvv2 = Gvv(1)**2 + Gvv(2)**2 + Gvv(3)**2
    IF( 0.5d0*Gvv2 <= ecutrho ) THEN 
      Ng = Ng + 1
    ENDIF 
  ENDDO
  ENDDO
  ENDDO

END FUNCTION 


END SUBROUTINE

