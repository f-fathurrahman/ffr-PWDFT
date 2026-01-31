REAL(8) function rfinp(rfmt1,rfir1,rfmt2,rfir2)
  USE modmain
  ! !INPUT/OUTPUT PARAMETERS:
  !   rfmt1 : first function in real spherical harmonics for all muffin-tins
  !           (in,real(npmtmax,natmtot))
  !   rfir1 : first real interstitial function in real-space (in,real(ngtot))
  !   rfmt2 : second function in real spherical harmonics for all muffin-tins
  !           (in,real(npmtmax,natmtot))
  !   rfir2 : second real interstitial function in real-space (in,real(ngtot))
  ! !DESCRIPTION:
  !   Calculates the inner product of two real functions over the entire unit
  !   cell. The input muffin-tin functions should have angular momentum cut-off
  !   {\tt lmaxo}. In the interstitial region, the integrand is multiplied with
  !   the characteristic function, $\tilde{\Theta}({\bf r})$, to remove the
  !   contribution from the muffin-tin. See routines {\tt rfmtinp} and
  !   {\tt gencfun}.
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(in) :: rfmt1(npmtmax,natmtot),rfir1(ngtot)
  REAL(8), intent(in) :: rfmt2(npmtmax,natmtot),rfir2(ngtot)
  ! local variables
  INTEGER :: is,ias,ir
  ! external functions
  REAL(8), external :: rfmtinp

  ! interstitial contribution
  rfinp=0.d0
  DO ir=1,ngtot
    rfinp=rfinp+rfir1(ir)*rfir2(ir)*cfunir(ir)
  ENDDO 
  rfinp=rfinp*omega/dble(ngtot)

  ! muffin-tin contribution
  DO ias=1,natmtot
    is=idxis(ias)
    rfinp=rfinp+rfmtinp(nrmt(is),nrmti(is),wrmt(:,is),rfmt1(:,ias),rfmt2(:,ias))
  ENDDO 
  RETURN 
END FUNCTION 
