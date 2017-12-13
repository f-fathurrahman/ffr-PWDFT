MODULE m_xc
  
  USE xc_f90_types_m
  USE xc_f90_lib_m
  IMPLICIT NONE 
  CHARACTER(10) :: XC_NAME = 'VWN'
  ! possible values: 'VWN' and 'PBE'

  TYPE(xc_f90_pointer_t) :: xc_func_x
  TYPE(xc_f90_pointer_t) :: xc_info_x

  TYPE(xc_f90_pointer_t) :: xc_func_c
  TYPE(xc_f90_pointer_t) :: xc_info_c

  LOGICAL :: USE_ARIAS_VWN = .FALSE.  ! only for debugging

END MODULE 


SUBROUTINE alloc_xc()
  USE m_xc
  IMPLICIT NONE 

  CALL xc_f90_func_init(xc_func_x, xc_info_x, 1, XC_UNPOLARIZED)
  CALL xc_f90_func_init(xc_func_c, xc_info_c, 1, XC_UNPOLARIZED)

END SUBROUTINE 

SUBROUTINE info_xc()
  USE m_xc
  IMPLICIT NONE 
END SUBROUTINE 

SUBROUTINE dealloc_xc()
  USE m_xc
  IMPLICIT NONE 
END SUBROUTINE 
