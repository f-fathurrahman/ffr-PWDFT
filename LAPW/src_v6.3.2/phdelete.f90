
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdelete
use modphonon
implicit none
! local variables
integer ios
character(256) fext
! construct the phonon file extension
call phfext(iqph,isph,iaph,ipph,fext)
! delete the eigenvector files
open(222,file=trim(scrpath)//'DEVECFV'//trim(fext),iostat=ios)
close(222,status='DELETE')
open(226,file=trim(scrpath)//'DEVECSV'//trim(fext),iostat=ios)
close(226,status='DELETE')
return
end subroutine

