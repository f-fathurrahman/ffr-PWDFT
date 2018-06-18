SUBROUTINE dealloc_PWGrid()

  USE m_PWGrid

  IF( allocated(Gv) ) DEALLOCATE( Gv )
  IF( allocated(Gv2) ) DEALLOCATE( Gv2 )
  IF( allocated(idx_g2r) ) DEALLOCATE( idx_g2r )
  IF( allocated(idx_gw2r) ) DEALLOCATE( idx_gw2r )
  IF( allocated(idx_gw2g) ) DEALLOCATE( idx_gw2g )

END SUBROUTINE 
