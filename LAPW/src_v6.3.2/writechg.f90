
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine writechg(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer is,ia,ias
! write charges
write(fnum,*)
write(fnum,'("Charges :")')
write(fnum,'(" core",T30,": ",G18.10)') chgcrtot
write(fnum,'(" valence",T30,": ",G18.10)') chgval
write(fnum,'(" interstitial",T30,": ",G18.10)') chgir
write(fnum,'(" muffin-tins (core leakage)")')
do is=1,nspecies
  write(fnum,'("  species : ",I4," (",A,")")') is,trim(spsymb(is))
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(fnum,'("   atom ",I4,T30,": ",G18.10," (",G18.10,")")') ia, &
     chgmt(ias),chgcrlk(ias)
  end do
end do
write(fnum,'(" total in muffin-tins",T30,": ",G18.10)') chgmttot
if (chgexs.ne.0.d0) then
  write(fnum,'(" excess",T30,": ",G18.10)') chgexs
end if
write(fnum,'(" total calculated charge",T30,": ",G18.10)') chgcalc
write(fnum,'(" total charge",T30,": ",G18.10)') chgtot
write(fnum,'(" error",T30,": ",G18.10)') abs(chgtot-chgcalc)
flush(fnum)
return
end subroutine

