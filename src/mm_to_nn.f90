FUNCTION mm_to_nn( mm, S ) RESULT(res)

  IMPLICIT NONE 
  INTEGER :: mm
  INTEGER :: S
  INTEGER :: res

  IF( mm > S/2 ) THEN 
    res = mm - S
    RETURN 
  ELSE
    res = mm
    RETURN 
  ENDIF 

END FUNCTION 

