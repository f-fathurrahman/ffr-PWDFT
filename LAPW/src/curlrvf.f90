SUBROUTINE curlrvf(rvfmt,rvfir,curlmt,curlir)
  USE modmain
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(in) :: rvfmt(npmtmax,natmtot,3),rvfir(ngtot,3)
  REAL(8), intent(out) :: curlmt(npmtmax,natmtot,3),curlir(ngtot,3)
  ! local variables
  INTEGER is,ias,np,i
  ! ALLOCATABLE arrays
  REAL(8), ALLOCATABLE :: grfmt(:,:,:,:),grfir(:,:,:)
  ALLOCATE(grfmt(npmtmax,natmtot,3,3),grfir(ngtot,3,3))

  ! compute the gradients
  DO i=1,3
    CALL gradrf(rvfmt(:,:,i),rvfir(:,i),grfmt(:,:,:,i),grfir(:,:,i))
  ENDDO 
  
  ! determine the muffin-tin and interstitial curl
  DO ias=1,natmtot
    is=idxis(ias)
    np=npmt(is)
    curlmt(1:np,ias,1)=grfmt(1:np,ias,2,3)-grfmt(1:np,ias,3,2)
    curlmt(1:np,ias,2)=grfmt(1:np,ias,3,1)-grfmt(1:np,ias,1,3)
    curlmt(1:np,ias,3)=grfmt(1:np,ias,1,2)-grfmt(1:np,ias,2,1)
  ENDDO 
  curlir(:,1)=grfir(:,2,3)-grfir(:,3,2)
  curlir(:,2)=grfir(:,3,1)-grfir(:,1,3)
  curlir(:,3)=grfir(:,1,2)-grfir(:,2,1)
  DEALLOCATE(grfmt,grfir)
  RETURN 
END SUBROUTINE 

