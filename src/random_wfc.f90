SUBROUTINE random_wfc( Nrows, Ncols, v )

  IMPLICIT NONE 
  INTEGER :: Nrows, Ncols
  COMPLEX(8) :: v(Nrows,Ncols)
  REAL(8) :: r1, r2
  INTEGER :: ic, ir

  DO ic = 1,Ncols
    DO ir = 1,Nrows
      CALL random_number( r1 )
      CALL random_number( r2 )
      v(ir,ic) = cmplx( r1, r2, kind=8 )
    ENDDO 
  ENDDO 

  CALL z_ortho_gram_schmidt( v, Nrows, Nrows, Ncols )

END SUBROUTINE 

