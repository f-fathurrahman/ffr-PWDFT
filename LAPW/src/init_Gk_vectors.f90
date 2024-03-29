SUBROUTINE init_Gk_vectors()
  USE m_atoms, ONLY: natmtot
  USE m_gvectors, ONLY: ngvec, vgc, ivg
  USE m_gkvectors, ONLY: ngkmax, gkmax, vgkc, sfacgk, ngk, igkig, vgkl, gkc
  USE m_kpoints, ONLY: vkc, nkpt, nkptnr, vkl
  USE m_spin, ONLY: vqcss, vqlss
  USE m_misc, ONLY: task
  USE m_spin, ONLY: spinsprl, nspnfv
  USE m_density_pot_xc, ONLY: xctype
  IMPLICIT NONE
  INTEGER :: ik, nppt, jspn
  REAL(8) :: vl(3), vc(3)

  !---------------------!
  !     G+k-vectors     !
  !---------------------!
  
  IF( (xctype(1) .lt. 0) .or. ( any(task==[5,10,205,300,600,620,630])) ) THEN 
    nppt = nkptnr
  ELSE 
    nppt = nkpt
  ENDIF 
  
  ! find the maximum number of G+k-vectors
  CALL findngkmax( nkpt, vkc, nspnfv, vqcss, ngvec, vgc, gkmax, ngkmax )
  
  ! allocate the G+k-vector arrays
  IF( allocated(ngk) ) DEALLOCATE(ngk)
  ALLOCATE( ngk(nspnfv,nppt) )
  !
  IF( allocated(igkig) ) DEALLOCATE(igkig)
  ALLOCATE( igkig(ngkmax,nspnfv,nppt) )
  !
  IF( allocated(vgkl)) DEALLOCATE(vgkl)
  ALLOCATE( vgkl(3,ngkmax,nspnfv,nppt) )
  !
  IF( allocated(vgkc)) DEALLOCATE(vgkc)
  ALLOCATE( vgkc(3,ngkmax,nspnfv,nppt) )
  !
  IF( allocated(gkc)) DEALLOCATE(gkc)
  ALLOCATE( gkc(ngkmax,nspnfv,nppt) )
  !
  IF( allocated(sfacgk)) DEALLOCATE(sfacgk)
  ALLOCATE( sfacgk(ngkmax, natmtot, nspnfv,nppt) )
  !
  DO ik=1,nppt
    DO jspn=1,nspnfv
      vl(:) = vkl(:,ik)
      vc(:) = vkc(:,ik)
      ! spin-spiral case
      IF( spinsprl ) THEN 
        IF(jspn == 1) THEN 
          vl(:) = vl(:) + 0.5d0*vqlss(:)
          vc(:) = vc(:) + 0.5d0*vqcss(:)
        ELSE 
          vl(:) = vl(:) - 0.5d0*vqlss(:)
          vc(:) = vc(:) - 0.5d0*vqcss(:)
        ENDIF 
      ENDIF 
      ! generate the G+k-vectors
      CALL gengkvec(ngvec,ivg,vgc,vl,vc,gkmax,ngkmax,ngk(jspn,ik), &
       igkig(:,jspn,ik),vgkl(:,:,jspn,ik),vgkc(:,:,jspn,ik),gkc(:,jspn,ik))
      ! generate structure factors for G+k-vectors
      CALL gensfacgp(ngk(jspn,ik),vgkc(:,:,jspn,ik),ngkmax,sfacgk(:,:,jspn,ik))
    ENDDO 
  ENDDO 

END SUBROUTINE 
