SUBROUTINE z_ortho_check( Npoints, Ncols, dVol, v )
  IMPLICIT NONE 
  INTEGER :: Npoints
  INTEGER :: Ncols
  REAL(8) :: dVol
  COMPLEX(8) :: v(Npoints,Ncols)
  !
  INTEGER :: ic
  COMPLEX(8) :: nrm
  !
  COMPLEX(8) :: zdotc

  WRITE(*,*)
  WRITE(*,*) 'Checking orthonormalization'
  WRITE(*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^'
  WRITE(*,*)
  WRITE(*,*) 'Norms:'
  DO ic = 1, Ncols
    nrm = zdotc( Npoints, v(:,ic), 1, v(:,ic),1 ) * dVol
    WRITE(*,'(1x,I8,2F18.10)') ic, nrm
  ENDDO
  
  WRITE(*,*)
  WRITE(*,*) 'Check wrt to vector 1'
  DO ic = 2,Ncols
    nrm = zdotc( Npoints, v(:,ic), 1, v(:,1),1 ) * dVol
    WRITE(*,'(1x,I8,2F18.10)') ic, nrm
  ENDDO

END SUBROUTINE 
