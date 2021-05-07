SUBROUTINE gradrf(rfmt,rfir,grfmt,grfir)
use modmain
IMPLICIT NONE 
! arguments
REAL(8), intent(in) :: rfmt(npmtmax,natmtot),rfir(ngtot)
REAL(8), intent(out) :: grfmt(npmtmax,natmtot,3),grfir(ngtot,3)
! local variables
INTEGER :: is,ias,ld,i
INTEGER :: ig,ifg
COMPLEX(8) :: z1
! ALLOCATABLE arrays
COMPLEX(8), ALLOCATABLE :: zfft1(:),zfft2(:)
! muffin-tin gradient
ld=npmtmax*natmtot

DO ias=1,natmtot
  is=idxis(ias)
  CALL gradrfmt(nrmt(is),nrmti(is),rlmt(:,1,is),rlmt(:,-1,is),rfmt(:,ias),ld, &
   grfmt(1,ias,1))
ENDDO 

! interstitial gradient
ALLOCATE(zfft1(ngtot),zfft2(ngtot))
zfft1(:)=rfir(:)
CALL zfftifc(3,ngridg,-1,zfft1)
DO i=1,3
  zfft2(:)=0.d0
  DO ig=1,ngvec
    ifg=igfft(ig)
    z1=zfft1(ifg)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(z1),dble(z1),8)
  ENDDO 
  CALL zfftifc(3,ngridg,1,zfft2)
  grfir(:,i)=dble(zfft2(:))
ENDDO 
DEALLOCATE(zfft1,zfft2)
RETURN 
END SUBROUTINE 

