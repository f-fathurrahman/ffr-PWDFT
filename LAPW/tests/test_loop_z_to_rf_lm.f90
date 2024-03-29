SUBROUTINE do_loop(lmax)
  IMPLICIT NONE 
  ! arguments
  INTEGER, INTENT(in) :: lmax
  ! local variables
  INTEGER :: l,m,lm1,lm2
  
  lm1 = 0
  DO l=0,lmax
    lm2 = lm1 + 2*(l+1)
    DO m=-l,-1
      lm1 = lm1 + 1
      lm2 = lm2 - 1
      write(*,'(1x,2I4,A,2I4)') l, m, ' : ', lm1, lm2
      !IF(mod(m,2) /= 0) THEN 
      !  rflm(lm1) = -c1*(aimag(zflm(lm1)) + aimag(zflm(lm2)))
      !ELSE 
      !  rflm(lm1) = c1*(aimag(zflm(lm2)) - aimag(zflm(lm1)))
      !ENDIF 
    ENDDO 
    lm1 = lm1+1
    lm2 = lm2-1
    !rflm(lm1) = dble(zflm(lm1))
    DO m = 1,l
      lm1 = lm1 + 1
      lm2 = lm2 - 1
      write(*,'(1x,2I4,A,2I4)') l, m, ' : ', lm1, lm2
      !IF( mod(m,2) /= 0) THEN 
      !  rflm(lm1) = c1*(dble(zflm(lm1)) - dble(zflm(lm2)))
      !ELSE 
      !  rflm(lm1) = c1*(dble(zflm(lm1)) + dble(zflm(lm2)))
      !ENDIF 
    ENDDO 
  ENDDO 
  RETURN 
END SUBROUTINE 


program main
  call do_loop(2)
end program

