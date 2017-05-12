SUBROUTINE calc_grad( Ncols, psi, grad )
  USE m_PWGrid, ONLY : Ngwx
  USE m_states, ONLY : Focc
  IMPLICIT NONE 
  INTEGER :: Ncols
  COMPLEX(8) :: psi(Ngwx,Ncols)
  COMPLEX(8) :: grad(Ngwx,Ncols)
  !
  COMPLEX(8), ALLOCATABLE :: Hpsi(:)
  INTEGER :: ic, icc
  !
  COMPLEX(8) :: zdotc

  ALLOCATE( Hpsi(Ngwx) )
  
  DO ic = 1,Ncols
    CALL op_H( 1, psi(:,ic), Hpsi )
    grad(:,ic) = Hpsi(:)
    DO icc = 1,Ncols
      grad(:,ic) = grad(:,ic) - zdotc( Ngwx, psi(:,icc),1, Hpsi(:),1 ) * psi(:,icc)
    ENDDO 
    ! multiply with Focc ???
    grad(:,ic) = Focc(ic)*grad(:,ic)
  ENDDO 

  DEALLOCATE( Hpsi )

END SUBROUTINE 
