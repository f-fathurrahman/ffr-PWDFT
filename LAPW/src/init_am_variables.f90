!-----------------------------
SUBROUTINE init_am_variables()
!-----------------------------
  USE m_constants, ONLY: zil, zilc, zi
  USE m_muffin_tins, ONLY: lmaxi, idxim, idxil, idxlm, lmmaxapw, lmmaxo, &
      lmaxo, lmmaxi, lmaxapw
  USE m_dos_optics_response, ONLY: lmaxdos
  IMPLICIT NONE 
  INTEGER :: l, m, lm

  WRITE(*,*) 'Setting up angular momentum variables'

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

END SUBROUTINE 

