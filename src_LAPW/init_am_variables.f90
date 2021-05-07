!-----------------------------
SUBROUTINE init_am_variables()
!-----------------------------
  USE modmain, ONLY: &
      lmaxi, zil, zilc, idxim, idxil, idxlm, lmmaxapw, lmmaxo, &
      lmaxo, zi, lmaxdos, lmmaxi, lmaxapw

  IMPLICIT NONE 
  INTEGER :: l, m, lm

  !------------------------------------!
  !     angular momentum variables     !
  !------------------------------------!
  IF( lmaxo > lmaxapw) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(init0): lmaxo > lmaxapw : ",2I8)') lmaxo,lmaxapw
    WRITE(*,*)
    STOP 
  ENDIF 
  lmaxi = min(lmaxi,lmaxo)
  lmmaxapw = (lmaxapw+1)**2
  lmmaxi = (lmaxi+1)**2
  lmmaxo = (lmaxo+1)**2
  ! check DOS lmax is within range
  lmaxdos = min(lmaxdos,lmaxo)

  ! index to (l,m) pairs
  IF( allocated(idxlm) ) DEALLOCATE(idxlm)
  ALLOCATE(idxlm(0:lmaxapw,-lmaxapw:lmaxapw))
  IF( allocated(idxil) ) DEALLOCATE(idxil)
  ALLOCATE( idxil(lmmaxapw) )
  IF( allocated(idxim) ) DEALLOCATE(idxim)
  ALLOCATE( idxim(lmmaxapw) )
  !
  lm = 0
  DO l = 0,lmaxapw
    DO m = -l,l
      lm = lm + 1
      idxlm(l,m) = lm
      idxil(lm) = l
      idxim(lm) = m
    ENDDO 
  ENDDO 

  ! array of i^l and (-i)^l values
  IF( allocated(zil) ) DEALLOCATE(zil)
  IF( allocated(zilc) ) DEALLOCATE(zilc)
  ALLOCATE( zil(0:lmaxapw), zilc(0:lmaxapw) )

  DO l = 0,lmaxapw
    zil(l) = zi**l
    zilc(l) = conjg(zil(l))
  ENDDO 
  
  ! write to VARIABLES.OUT
  WRITE(*,*) 'lmaxapw = ', lmaxapw
  WRITE(*,*) 'lmaxi = ', lmaxi
  WRITE(*,*) 'lmaxo = ', lmaxo

END SUBROUTINE 

