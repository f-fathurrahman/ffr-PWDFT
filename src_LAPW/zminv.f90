
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

SUBROUTINE zminv(n,a)
IMPLICIT NONE 
! arguments
INTEGER, intent(in) :: n
COMPLEX(8), intent(inout) :: a(n,n)
! local variables
INTEGER :: info
! automatic arrays
INTEGER :: ipiv(n)
COMPLEX(8) work(n)

CALL zgetrf(n,n,a,n,ipiv,info)
IF(info.ne.0) THEN 
  WRITE(*,*)
  WRITE(*,'("Error(zminv): unable to invert matrix")')
  WRITE(*,'(" ZGETRF RETURN ed INFO = ",I8)') info
  WRITE(*,*)
  stop
ENDIF 
CALL zgetri(n,a,n,ipiv,work,n,info)
IF(info.ne.0) THEN 
  WRITE(*,*)
  WRITE(*,'("Error(zminv): unable to invert matrix")')
  WRITE(*,'(" ZGETRI RETURN ed INFO = ",I8)') info
  WRITE(*,*)
  stop
ENDIF 

RETURN 
END SUBROUTINE 

