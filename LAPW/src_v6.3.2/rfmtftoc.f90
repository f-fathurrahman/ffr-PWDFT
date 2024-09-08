
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfmtftoc(nrc,nrci,rfmt,rfcmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nrc,nrci
real(8), intent(in) :: rfmt(*)
real(8), intent(out) :: rfcmt(*)
! local variables
integer irc,i,j,n
i=1
j=1
n=lmmaxi*lradstp
do irc=1,nrci
  call dcopy(lmmaxi,rfmt(i),1,rfcmt(j),1)
  i=i+n
  j=j+lmmaxi
end do
i=i+(lradstp-1)*(lmmaxo-lmmaxi)
n=lmmaxo*lradstp
do irc=nrci+1,nrc
  call dcopy(lmmaxo,rfmt(i),1,rfcmt(j),1)
  i=i+n
  j=j+lmmaxo
end do
return
end subroutine

