!---------------------------------
SUBROUTINE init_vector_field_E_A()
!---------------------------------
  USE modmain, ONLY: ainv, tshift, tevecsv, tefield, tafield, &
                     efieldc, efieldl, epslat, afieldc, afieldl
  IMPLICIT NONE 

  !-------------------------------!
  !     vector fields E and A     !
  !-------------------------------!
  tefield=.false.
  IF( sum(abs(efieldc(:))) > epslat) THEN 
    ! no shift of the atomic positions
    tshift = .false.
    ! electric field vector in lattice coordinates
    call r3mv( ainv, efieldc, efieldl )
    tefield = .true.
  ENDIF 
  
  tafield = .false.
  IF( sum(abs(afieldc(:))) > epslat ) THEN 
    tafield = .true.
    ! A-field in lattice coordinates
    call r3mv(ainv, afieldc, afieldl)
    ! vector potential added in second-variational step
    tevecsv = .true.
  ENDIF 

END SUBROUTINE 
