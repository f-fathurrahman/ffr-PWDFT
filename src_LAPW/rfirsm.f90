SUBROUTINE rfirsm(m,rfir)
  USE m_gvectors, ONLY: ngtot, ngridg, gc, gmaxvr, igfft
  IMPLICIT NONE 
  ! arguments
  INTEGER, intent(in) :: m
  REAL(8), intent(inout) :: rfir(ngtot)
  ! local variables
  INTEGER :: ig,ifg
  REAL(8) :: t0,t1,t2
  ! ALLOCATABLE arrays
  COMPLEX(8), ALLOCATABLE :: zfft(:)
  IF(m.le.0) RETURN 
  ALLOCATE(zfft(ngtot))
  zfft(:)=rfir(:)
  CALL zfftifc(3,ngridg,-1,zfft)
  t0=dble(2*m)
  t1=1.d0/gmaxvr
  DO ig=1,ngtot
    ifg=igfft(ig)
    t2=t1*gc(ig)
    t2=exp(-t0*t2**4)
    zfft(ifg)=t2*zfft(ifg)
  ENDDO 
  CALL zfftifc(3,ngridg,1,zfft)
  rfir(:)=dble(zfft(:))
  DEALLOCATE(zfft)
  RETURN 
END SUBROUTINE 

