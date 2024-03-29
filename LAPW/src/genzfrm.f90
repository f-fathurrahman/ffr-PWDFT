SUBROUTINE genzfrm(wfmt11,wfmt12,wfir11,wfir12,wfmt21,wfmt22,wfir21,wfir22, &
 zrhomt,zrhoir,zmagmt,zmagir)
  use modmain
  IMPLICIT NONE 
  ! arguments
  COMPLEX(8), intent(in) ::  wfmt11(npcmtmax,natmtot),wfmt12(npcmtmax,natmtot)
  COMPLEX(8), intent(in) ::  wfir11(ngtot),wfir12(ngtot)
  COMPLEX(8), intent(in) ::  wfmt21(npcmtmax,natmtot),wfmt22(npcmtmax,natmtot)
  COMPLEX(8), intent(in) ::  wfir21(ngtot),wfir22(ngtot)
  COMPLEX(8), intent(out) :: zrhomt(npcmtmax,natmtot),zrhoir(ngtot)
  COMPLEX(8), intent(out) :: zmagmt(npcmtmax,natmtot,ndmag),zmagir(ngtot,ndmag)
  ! local variables
  INTEGER ld,is,ias,nthd
  !-------------------------!
  !     muffin-tin part     !
  !-------------------------!
  ld=npcmtmax*natmtot
  DO ias=1,natmtot
    is=idxis(ias)
    CALL genzrm(npcmt(is),wfmt11(:,ias),wfmt12(:,ias),wfmt21(:,ias), &
     wfmt22(:,ias),zrhomt(:,ias),ld,zmagmt(:,ias,1))
  ENDDO 
  !---------------------------!
  !     interstitial part     !
  !---------------------------!
  CALL genzrm(ngtot,wfir11,wfir12,wfir21,wfir22,zrhoir,ngtot,zmagir)
  RETURN 
END SUBROUTINE 

