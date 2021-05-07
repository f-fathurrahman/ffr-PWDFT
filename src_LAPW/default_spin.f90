SUBROUTINE default_spin()
  USE m_spin, ONLY: ssdph, spinsprl, spinpol, spinorb, socscf, rndbfcmt, reducebf, npmae0, &
              fsmtype, cmagz, vqlss, mommtfix, momfix, bfieldc0, bfcmt0
  IMPLICIT NONE 

  spinsprl=.false.
  ssdph=.true.
  vqlss(:)=0.d0
  spinpol=.false.
  spinorb=.false.
  rndbfcmt=0.d0
  bfcmt0(:,:,:)=0.d0
  fsmtype=0
  momfix(:)=0.d0
  mommtfix(:,:,:)=1.d6
  cmagz=.false.
  npmae0=-1
  socscf=1.d0
  reducebf=1.d0
  bfieldc0(:)=0.d0

END SUBROUTINE 