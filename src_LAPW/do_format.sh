#!/bin/bash
sed -i "s/do /DO /g" "$1"
sed -i "s/end do/ENDDO /g" "$1"
sed -i "s/else if (/ELSEIF(/g" "$1"
sed -i "s/if (/IF(/g" "$1"
sed -i "s/end if/ENDIF /g" "$1"
sed -i "s/then/THEN /g" "$1"
sed -i "s/end subroutine/END SUBROUTINE /g" "$1"
sed -i "s/subroutine /SUBROUTINE /g" "$1"
sed -i "s/return/RETURN /g" "$1"
sed -i "s/integer /INTEGER /g" "$1"
sed -i "s/integer, /INTEGER, /g" "$1"
sed -i "s/logical /LOGICAL /g" "$1"
sed -i "s/logical, /LOGICAL, /g" "$1"
sed -i "s/real(8)/REAL(8)/g" "$1"
sed -i "s/complex(8)/COMPLEX(8)/g" "$1"
sed -i "s/implicit none/IMPLICIT NONE /g" "$1"
sed -i "s/call /CALL /g" "$1"
sed -i "s/allocatable /ALLOCATABLE /g" "$1"
sed -i "s/deallocate(/DEALLOCATE(/g" "$1"
sed -i "s/allocate(/ALLOCATE(/g" "$1"
sed -i "s/write(/WRITE(/g" "$1"