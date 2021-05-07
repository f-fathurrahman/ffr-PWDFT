COMPLEX(8) function zfinp(zfmt1,zfir1,zfmt2,zfir2)
  ! !USES:
  use modmain
! !INPUT/OUTPUT PARAMETERS:
!   zfmt1 : first complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(npcmtmax,natmtot))
!   zfir1 : first complex interstitial function in real-space
!           (in,complex(ngtc))
!   zfmt2 : second complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(npcmtmax,natmtot))
!   zfir2 : second complex interstitial function in real-space
!           (in,complex(ngtc))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions over the entire unit
!   cell. The muffin-tin functions should be stored on the coarse radial grid.
!   In the interstitial region, the integrand is multiplied with the
!   characteristic function to remove the contribution from the muffin-tin. See
!   routines {\tt zfmtinp} and {\tt gencfun}.
  IMPLICIT NONE 
  ! arguments
  COMPLEX(8), intent(in) :: zfmt1(npcmtmax,natmtot),zfir1(ngtc)
  COMPLEX(8), intent(in) :: zfmt2(npcmtmax,natmtot),zfir2(ngtc)
  ! local variables
  INTEGER :: is,ias,ir,nthd
  ! external functions
  COMPLEX(8) zfmtinp
  external zfmtinp

  ! interstitial contribution
  zfinp=0.d0
  DO ir=1,ngtc
    zfinp=zfinp+cfrc(ir)*conjg(zfir1(ir))*zfir2(ir)
  ENDDO 
  zfinp=zfinp*(omega/dble(ngtc))
  ! muffin-tin contribution
  DO ias=1,natmtot
    is=idxis(ias)
    zfinp=zfinp+zfmtinp(nrcmt(is),nrcmti(is),wrcmt(:,is),zfmt1(:,ias), &
     zfmt2(:,ias))
  ENDDO 
  RETURN 
END FUNCTION 