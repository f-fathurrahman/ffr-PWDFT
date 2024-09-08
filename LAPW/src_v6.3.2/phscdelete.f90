
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phscdelete
use modmain
implicit none
! local variables
integer ios
! delete the eigenvector files
call delevec
! delete the eigenvalue files
open(120,file='EVALFV'//trim(filext),iostat=ios)
close(120,status='DELETE')
open(124,file='EVALSV'//trim(filext),iostat=ios)
close(124,status='DELETE')
! delete the occupancy file
open(130,file='OCCSV'//trim(filext),iostat=ios)
close(130,status='DELETE')
return
end subroutine

