SUBROUTINE default_plotting()
  USE m_plotting, ONLY: vclp3d, vclp2d, nvp1d, npp1d, np3d, np2d, vvlp1d
  IMPLICIT NONE

  nvp1d=2
  IF(allocated(vvlp1d)) DEALLOCATE(vvlp1d)
  ALLOCATE(vvlp1d(3,nvp1d))
  vvlp1d(:,1)=0.d0
  vvlp1d(:,2)=1.d0
  npp1d=200
  vclp2d(:,:)=0.d0
  vclp2d(1,1)=1.d0
  vclp2d(2,2)=1.d0
  np2d(:)=40
  vclp3d(:,:)=0.d0
  vclp3d(1,1)=1.d0
  vclp3d(2,2)=1.d0
  vclp3d(3,3)=1.d0
  np3d(:)=20 
END SUBROUTINE