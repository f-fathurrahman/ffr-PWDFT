PROGRAM test_ortho

  IMPLICIT NONE 
  INTEGER :: Ngwx, Nstates, Ngw
  COMPLEX(8), ALLOCATABLE :: v(:,:)
  INTEGER :: ist, igw
  REAL(8) :: r1, r2

  Ngwx = 1000
  Ngw = 980
  Nstates = 4
  ALLOCATE( v(Ngwx,Nstates) )
  DO ist = 1,Nstates
    DO igw = 1,Ngw
      CALL random_number(r1)
      CALL random_number(r2)
      v(igw,ist) = cmplx( r1, r2, kind=8 )
    ENDDO 
  ENDDO 

  !CALL z_ortho_gram_schmidt( v, Ngwx, Ngwx, Nstates )
  CALL z_ortho_qr( v, Ngwx, Ngw, Nstates )

  CALL z_ortho_check( Ngw, Nstates, 1.d0, v(1:Ngw,:) )

  DEALLOCATE( v )

END PROGRAM 

