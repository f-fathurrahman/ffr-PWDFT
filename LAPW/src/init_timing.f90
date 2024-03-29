SUBROUTINE init_timing()
  USE m_timing
  IMPLICIT NONE 
  !-------------------------------!
  !     zero timing variables     !
  !-------------------------------!
  timeinit=0.d0
  timemat=0.d0
  timefv=0.d0
  timesv=0.d0
  timerho=0.d0
  timepot=0.d0
  timefor=0.d0
END SUBROUTINE 
