SUBROUTINE moment()
use modmain
! !DESCRIPTION:
!   Computes the muffin-tin, interstitial and total moments by integrating the
!   magnetisation.
IMPLICIT NONE 
! local variables
INTEGER idm,is,ias,nr,nri
REAL(8) t1
! automatic arrays
REAL(8) fr(nrmtmax)
IF(.not.spinpol) THEN 
  mommt(:,:)=0.d0
  mommttot(:)=0.d0
  momir(:)=0.d0
  momtot(:)=0.d0
  RETURN 
ENDIF 
! find the muffin-tin moments
mommttot(:)=0.d0
DO idm=1,ndmag
  DO ias=1,natmtot
    is=idxis(ias)
    nr=nrmt(is)
    nri=nrmti(is)
! extract the l=m=0 component from the muffin-tin magnetisation
    CALL rfmtlm(1,nr,nri,magmt(:,ias,idm),fr)
! integrate to the muffin-tin radius
    t1=dot_product(wrmt(1:nr,is),fr(1:nr))
    mommt(idm,ias)=fourpi*y00*t1
    mommttot(idm)=mommttot(idm)+mommt(idm,ias)
  ENDDO 
ENDDO 
! find the interstitial moments
DO idm=1,ndmag
  t1=dot_product(magir(:,idm),cfunir(:))
  momir(idm)=t1*omega/dble(ngtot)
ENDDO 
momtot(:)=mommttot(:)+momir(:)
! total moment magnitude
IF(ncmag) THEN 
  momtotm=sqrt(momtot(1)**2+momtot(2)**2+momtot(3)**2)
ELSE 
  momtotm=abs(momtot(1))
ENDIF 
! write total moment magnitude to test file
WRITE(*,*) 'total moment magnitude', momtotm
RETURN 
END SUBROUTINE 
