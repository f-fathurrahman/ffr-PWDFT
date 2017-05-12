PROGRAM test_ortho

  IMPLICIT NONE 
  INTEGER :: Ngwx, Nstates
  COMPLEX(8), ALLOCATABLE :: v(:,:)
  INTEGER :: ist, igw
  REAL(8) :: r1, r2

  Ngwx = 100
  Nstates = 4
  ALLOCATE( v(Ngwx,Nstates) )
  DO ist = 1,Nstates
    DO igw = 1,Ngwx
      CALL random_number(r1)
      CALL random_number(r2)
      v(igw,ist) = cmplx( r1, r2, kind=8 )
    ENDDO 
  ENDDO 

  CALL z_ortho_gram_schmidt( v, Ngwx, Ngwx, Nstates )

  CALL z_ortho_check( Ngwx, Nstates, 1.d0, v )

  DEALLOCATE( v )

END PROGRAM 

