SUBROUTINE writeengy(fnum)
use modmain
use moddftu
IMPLICIT NONE 
! arguments
INTEGER, intent(in) :: fnum
WRITE(fnum,*)
WRITE(fnum,'("Energies :")')
WRITE(fnum,'(" Fermi",T30,": ",F22.12)') efermi
WRITE(fnum,'(" sum of eigenvalues",T30,": ",F22.12)') evalsum
WRITE(fnum,'(" electron kinetic",T30,": ",F22.12)') engykn
WRITE(fnum,'(" core electron kinetic",T30,": ",F22.12)') engykncr
WRITE(fnum,'(" Coulomb",T30,": ",F22.12)') engycl
WRITE(fnum,'(" Coulomb potential",T30,": ",F22.12)') engyvcl
WRITE(fnum,'(" nuclear-nuclear",T30,": ",F22.12)') engynn
WRITE(fnum,'(" electron-nuclear",T30,": ",F22.12)') engyen
WRITE(fnum,'(" Hartree",T30,": ",F22.12)') engyhar
WRITE(fnum,'(" Madelung",T30,": ",F22.12)') engymad
WRITE(fnum,'(" xc potential",T30,": ",F22.12)') engyvxc
IF(spinpol) THEN 
  WRITE(fnum,'(" xc effective B-field",T30,": ",F22.12)') engybxc
  WRITE(fnum,'(" external B-field",T30,": ",F22.12)') engybext
ENDIF 

WRITE(fnum,'(" exchange",T30,": ",F22.12)') engyx
WRITE(fnum,'(" correlation",T30,": ",F22.12)') engyc
WRITE(fnum,'(" exchange-correlation",T30,": ",F22.12)') engyx + engyc
if (dftu.ne.0) then
  write(fnum,'(" DFT+U",T30,": ",F22.12)') engydu
end if
IF(stype.eq.3) THEN 
  WRITE(fnum,'(" electron entropic",T30,": ",F22.12)') engyts
ENDIF 
WRITE(fnum,'(" total energy",T30,": ",F22.12)') engytot
IF(spinpol) THEN 
  WRITE(fnum,'(" (external B-field energy excluded from total)")')
ENDIF 
flush(fnum)
RETURN 
END SUBROUTINE 

