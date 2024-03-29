COMPLEX(8) function rzfinp(rfmt,rfir,zfmt,zfir)
  use modmain
  IMPLICIT NONE 
  ! arguments
  REAL(8), intent(in) :: rfmt(npcmtmax,natmtot),rfir(ngtot)
  COMPLEX(8), intent(in) :: zfmt(npcmtmax,natmtot),zfir(ngtot)
  ! local variables
  INTEGER :: is,ias,ir
  ! external functions
  COMPLEX(8) rzfmtinp
  external rzfmtinp

  ! interstitial contribution
  rzfinp=0.d0
  DO ir=1,ngtot
    rzfinp=rzfinp+(cfunir(ir)*rfir(ir))*zfir(ir)
  ENDDO 
  rzfinp=rzfinp*(omega/dble(ngtot))

  ! muffin-tin contribution
  DO ias=1,natmtot
    is=idxis(ias)
    rzfinp=rzfinp+rzfmtinp(nrcmt(is),nrcmti(is),wrcmt(:,is),rfmt(:,ias), &
     zfmt(:,ias))
  ENDDO 
  RETURN 
END FUNCTION 

