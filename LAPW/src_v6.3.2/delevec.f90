
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: delevec
! !INTERFACE:
subroutine delevec
! !USES:
use modmain
! !DESCRIPTION:
!   Deletes the first- and second-variational eigenvector files {\tt EVECFV.OUT}
!   and {\tt EVECSV.OUT}.
!
! !REVISION HISTORY:
!   Created May 2007 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ios
! delete the first-variational eigenvector file
open(122,file=trim(scrpath)//'EVECFV'//trim(filext),iostat=ios)
close(122,status='DELETE')
! delete the second-variational eigenvector file
open(126,file=trim(scrpath)//'EVECSV'//trim(filext),iostat=ios)
close(126,status='DELETE')
return
end subroutine
!EOC

