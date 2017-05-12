SUBROUTINE init_rgrid()
  USE m_cell, ONLY : LatVecs
  USE m_realspace, ONLY : Ns, rgrid, Npoints
  IMPLICIT NONE 
  INTEGER :: ip, i, j, k

  ALLOCATE( rgrid(3,Npoints) )

  ip = 0
  DO k = 0, Ns(3)-1
  DO j = 0, Ns(2)-1
  DO i = 0, Ns(1)-1
    ip = ip + 1
    rgrid(1,ip) = LatVecs(1,1)*i/Ns(1) + LatVecs(2,1)*j/Ns(2) + LatVecs(3,1)*k/Ns(3)
    rgrid(2,ip) = LatVecs(1,2)*i/Ns(1) + LatVecs(2,2)*j/Ns(2) + LatVecs(3,2)*k/Ns(3)
    rgrid(3,ip) = LatVecs(1,3)*i/Ns(1) + LatVecs(2,3)*j/Ns(2) + LatVecs(3,3)*k/Ns(3)
  ENDDO 
  ENDDO 
  ENDDO 


END SUBROUTINE 

