subroutine nuclei
use modmain, only: nrcmt, nrcnucl, rnucl, nrmt, nrnucl, rsp, rcmt, &
                   volnucl, pi, spzn, nspecies
implicit none
! local variables
integer is,ir,irc
! external functions
real(8) radnucl
external radnucl
do is=1,nspecies
! approximate nuclear radius
  rnucl(is)=radnucl(spzn(is))
! nuclear volume
  volnucl(is)=(4.d0/3.d0)*pi*rnucl(is)**3
! number of radial mesh points to nuclear radius
  nrnucl(is)=1
  do ir=1,nrmt(is)
    if (rsp(ir,is).gt.rnucl(is)) then
      nrnucl(is)=ir
      exit
    end if
  end do
! number of coarse radial mesh points to nuclear radius
  nrcnucl(is)=1
  do irc=1,nrcmt(is)
    if (rcmt(irc,is).gt.rnucl(is)) then
      nrcnucl(is)=irc
      exit
    end if
  end do
end do
return
end subroutine

