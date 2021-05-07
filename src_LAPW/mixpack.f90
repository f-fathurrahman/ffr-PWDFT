SUBROUTINE mixpack(tpack,n,v)
  USE modmain
! !INPUT/OUTPUT PARAMETERS:
!   tpack : .true. for packing, .false. for unpacking (in,logical)
!   n     : total number of real values stored (out,integer)
!   v     : packed potential (inout,real(*))
! !DESCRIPTION:
!   Packs/unpacks the muffin-tin and interstitial parts of the Kohn-Sham
!   potential and magnetic field into/from the single array {\tt v}. This array
!   can THEN  be passed directly to the mixing routine. See routine {\tt rfpack}.
  IMPLICIT NONE 
  ! arguments
  LOGICAL, intent(in) :: tpack
  INTEGER, intent(out) :: n
  REAL(8), intent(inout) :: v(*)
  ! local variables
  INTEGER :: idm,ispn,jspn
  INTEGER :: ias,lm1,lm2
  n=0
  ! pack the Kohn-Sham potential and magnetic field
  CALL rfpack(tpack,n,npmt,npmtmax,vsmt,vsir,v)
  DO idm=1,ndmag
    CALL rfpack(tpack,n,npcmt,npcmtmax,bsmt(:,:,idm),bsir(:,idm),v)
  ENDDO 

  RETURN 
END SUBROUTINE 
