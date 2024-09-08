
! Copyright (C) 2016 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnulr(ik0,evecu)
use modmain
use modulr
implicit none
! arguments
integer, intent(in) :: ik0
complex(8), intent(out) :: evecu(nstulr,nstulr)
! generate the ultra long-range Hamiltonian
call genhmlu(ik0,evecu)
! find the eigenvalues and vectors
call eveqnz(nstulr,nstulr,evecu,evalu(:,ik0))
return
end subroutine

