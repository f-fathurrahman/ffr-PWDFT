SUBROUTINE dealloc_PsPot()

  USE m_PsPot
  IMPLICIT NONE 

  IF( allocated(PsPot_FilePath) ) DEALLOCATE(PsPot_FilePath)
  IF( allocated(Ps_HGH_Params) ) DEALLOCATE(Ps_HGH_Params)

  IF( allocated(w_NL) ) DEALLOCATE(w_NL)
  IF( allocated(betaNL) ) DEALLOCATE(betaNL)

  IF( allocated(SpecNprojTot) ) DEALLOCATE(SpecNprojTot)

  IF( allocated(PsPot_lll) ) DEALLOCATE(PsPot_lll)
  IF( allocated(PsPot_ipr) ) DEALLOCATE(PsPot_ipr)

END SUBROUTINE 

