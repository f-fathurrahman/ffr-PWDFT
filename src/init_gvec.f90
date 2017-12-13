SUBROUTINE init_gvec()
  USE m_constants, ONLY : PI
  USE m_realspace, ONLY : Ns, Npoints
  USE m_cell, ONLY : RecVecs
  USE m_PWGrid, ONLY : Gv, Gv2, Ng
  IMPLICIT NONE
  !
  INTEGER :: i, j, k, ig, ii, jj, kk
  ! Function
  INTEGER :: mm_to_nn

  Ng = Npoints
  ALLOCATE( Gv(3,Ng) )
  ALLOCATE( Gv2(Ng) )

  ig = 0
  DO k = 0, Ns(3)-1
  DO j = 0, Ns(2)-1
  DO i = 0, Ns(1)-1
    !
    ig = ig + 1
    !
    ii = mm_to_nn( i, Ns(1) )
    jj = mm_to_nn( j, Ns(2) )
    kk = mm_to_nn( k, Ns(3) )
    !
    Gv(1,ig) = RecVecs(1,1)*ii + RecVecs(2,1)*jj + RecVecs(3,1)*kk
    Gv(2,ig) = RecVecs(1,2)*ii + RecVecs(2,2)*jj + RecVecs(3,2)*kk
    Gv(3,ig) = RecVecs(1,3)*ii + RecVecs(2,3)*jj + RecVecs(3,3)*kk
    !
    Gv2(ig) = Gv(1,ig)**2 + Gv(2,ig)**2 + Gv(3,ig)**2
  ENDDO
  ENDDO
  ENDDO

END SUBROUTINE

