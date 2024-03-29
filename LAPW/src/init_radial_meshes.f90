!------------------------------
SUBROUTINE init_radial_meshes()
!------------------------------
  USE m_atoms, ONLY: nspecies
  USE m_muffin_tins, ONLY: nrmtmax, nrcmtmax, lradstp, nrcmt, nrcmti, &
                     npcmti, npcmtmax, npcmt, nrmti, npmt, nrmti, npmti, &
                     nrmt, npcmt, nrcmti, lmmaxi, lmmaxo, npmtmax
  IMPLICIT NONE 
  INTEGER :: is

  WRITE(*,*) 'Setting up radial meshes'

  nrmtmax = 1
  nrcmtmax = 1
  DO is = 1,nspecies
    ! make the muffin-tin mesh commensurate with lradstp
    nrmt(is) = nrmt(is) - mod(nrmt(is)-1, lradstp)
    nrmtmax = max(nrmtmax,nrmt(is))
    ! number of coarse radial mesh points
    nrcmt(is) = (nrmt(is)-1)/lradstp + 1
    nrcmtmax = max(nrcmtmax,nrcmt(is))
  ENDDO 

  ! set up atomic and muffin-tin radial meshes
  CALL genrmesh()

  ! number of points in packed muffin-tins
  npmtmax = 1
  npcmtmax = 1
  DO is=1,nspecies
    npmti(is) = lmmaxi*nrmti(is)
    npmt(is) = npmti(is) + lmmaxo*(nrmt(is)-nrmti(is))
    npmtmax = max(npmtmax,npmt(is))
    npcmti(is) = lmmaxi*nrcmti(is)
    npcmt(is) = npcmti(is) + lmmaxo*(nrcmt(is)-nrcmti(is))
    npcmtmax = max(npcmtmax, npcmt(is))
  ENDDO 

END SUBROUTINE 
