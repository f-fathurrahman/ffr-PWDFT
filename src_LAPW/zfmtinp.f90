COMPLEX(8) FUNCTION zfmtinp(nr,nri,wr,zfmt1,zfmt2)
  USE m_muffin_tins, ONLY: lmmaxi, lmmaxo, lmaxi
  ! !INPUT/OUTPUT PARAMETERS:
  !   nr    : number of radial mesh points (in,integer)
  !   nri   : number of points on the inner part of the muffin-tin (in,integer)
  !   wr    : weights for integration on radial mesh (in,real(nr))
  !   zfmt1 : first complex muffin-tin function in spherical harmonics
  !           (in,complex(*))
  !   zfmt2 : second complex muffin-tin function (in,complex(*))
  ! !DESCRIPTION:
  !   Calculates the inner product of two complex fuctions in the muffin-tin. In
  !   other words, given two complex functions of the form
  !   $$ f({\bf r})=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}f_{lm}(r)Y_{lm}
  !    (\hat{\bf r}), $$
  !   the function RETURN s
  !   $$ I=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}\int f_{lm}^{1*}(r)
  !    f_{lm}^2(r)r^2\,dr\;. $$
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: nr,nri
  REAL(8), intent(in) :: wr(nr)
  COMPLEX(8), intent(in) :: zfmt1(*),zfmt2(*)
  ! local variables
  INTEGER :: ir,i
  COMPLEX(8) :: z1
  ! automatic arrays
  REAL(8) :: fr1(nr),fr2(nr)
  ! external functions
  COMPLEX(8) zdotc
  external zdotc

  ! compute the dot-products for each radial point
  i=1
  IF(lmaxi.eq.1) THEN 
    DO ir=1,nri
      z1=conjg(zfmt1(i))*zfmt2(i) &
        +conjg(zfmt1(i+1))*zfmt2(i+1) &
        +conjg(zfmt1(i+2))*zfmt2(i+2) &
        +conjg(zfmt1(i+3))*zfmt2(i+3)
      fr1(ir)=dble(z1); fr2(ir)=aimag(z1)
      i=i+4
    ENDDO 
  ELSE 
    DO ir=1,nri
      z1=zdotc(lmmaxi,zfmt1(i),1,zfmt2(i),1)
      fr1(ir)=dble(z1); fr2(ir)=aimag(z1)
      i=i+lmmaxi
    ENDDO 
  ENDIF
  ! 
  DO ir=nri+1,nr
    z1=zdotc(lmmaxo,zfmt1(i),1,zfmt2(i),1)
    fr1(ir)=dble(z1); fr2(ir)=aimag(z1)
    i=i+lmmaxo
  ENDDO 
  ! integrate over r
  zfmtinp=cmplx(dot_product(wr(:),fr1(:)),dot_product(wr(:),fr2(:)),8)
  RETURN 
END FUNCTION 